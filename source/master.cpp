#include <fstream>
#include <cmath>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <math_constants.hpp>
#include <progress_meter.hpp>
#include <wigner_3j.hpp>
#include <master.hpp>
#include <matrix_impl.hpp>

#include <healpix_base.h>
#include <alm.h>
#include <alm_healpix_tools.h>
#include <xcomplex.h>
#include <fitshandle.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>


double lFactor(int l)
{
    if(l == 0)
        return 1.0 / (2 * Math::pi);
    else
        return double(l) * (l + 1) / (2 * Math::pi);
}

Master::Master(const Healpix_Map<double>& mask, const char* couplingKernelFileName, const std::vector<double>& beam, const std::vector<int>* bins, int lMax) : beam_(beam), couplingKernelFileName_(couplingKernelFileName), map_(NULL)
{
    output_screen("Copying the mask..." << std::endl);
    mask_.SetNside(mask.Nside(), RING);
    
    mask_.Import(mask);
    output_screen("OK" << std::endl);
    
    construct(bins, lMax);
}

Master::Master(const char* maskName, const char* couplingKernelFileName, const std::vector<double>& beam, const std::vector<int>* bins, int lMax) : beam_(beam), couplingKernelFileName_(couplingKernelFileName), map_(NULL)
{
    output_screen("Reading the mask..." << std::endl);
    read_Healpix_map_from_fits(std::string(maskName), mask_);
    
    if(mask_.Scheme() == NEST)
        mask_.swap_scheme();
    
    check(mask_.Scheme() == RING, "");
    
    output_screen("OK" << std::endl);
    
    construct(bins, lMax);
}

double absSq(const xcomplex<double>& x)
{
    const double re = x.real();
    const double im = x.imag();
    return re * re + im * im;
}

double ps(const Alm<xcomplex<double> >& alm, int l)
{
    check(l >= 0 && l <= alm.Lmax(), "invalid l");
    double res = absSq(alm(l, 0));
    for(int m = 1; m <= l; ++m)
        res += 2 * absSq(alm(l, m));
    
    res /= (2 * l + 1);
    
    return res;
}

void
Master::construct(const std::vector<int>* bins, int lMax)
{
    if(bins)
    {
        check(lMax == 0, "cannot give both a bin vector and lMax");
        lMax_ = 0;

        bins_ = *bins;
        check(bins_.size() > 1, "");
        check(bins_[0] >= 0, "");
#ifdef CHECKS_ON
        for(int i = 1; i < bins_.size(); ++i)
        {
            check(bins_[i] > bins_[i - 1], "");
        }
#endif
        lMax_ = bins_[bins_.size() - 1] - 1;
        check(lMax_ <= 10000, "lMax is too big");
        check(beam_.size() >= lMax_ + 1, "incomplete beam function");
    }
    else
    {
        check(lMax > 0, "invalid lMax = " << lMax);
        lMax_ = lMax;
    }

    output_screen("Calculating f_sky..." << std::endl);
    ProgressMeter meter(mask_.Npix());
    double fSky = 0;
    for(long i = 0; i < mask_.Npix(); ++i)
    {
        fSky += mask_[i];
        meter.advance();
    }
    fSky /= mask_.Npix();
    output_screen("OK" << std::endl);
    output_screen("f_sky = " << fSky << std::endl);
    
    Alm<xcomplex<double> > almMask(lMax_, lMax_);
    arr<double> weight(2 * mask_.Nside(), 1);
    
    output_screen("Calculating alm..." << std::endl);
    map2alm(mask_, almMask, weight);
    output_screen("OK" << std::endl);
    
    
    //pseudo power spectra
    check(w_.empty(), "");
    
    for(int l = 0; l <= lMax_; ++l)
        w_.push_back(ps(almMask, l));

    calculateCoupling();
    calculateK();
}

void
Master::calculate()
{
    check(map_, "");
    check(map_->Nside() == mask_.Nside(), "the map and the mask must have the same NSide");
    
    maskedMap_.SetNside(map_->Nside(), RING);
    
    for(long i = 0; i < map_->Npix(); ++i)
        maskedMap_[i] = (*map_)[i] * mask_[i];
    
    Alm<xcomplex<double> > alm(lMax_, lMax_);
    arr<double> weight(2 * map_->Nside(), 1);

    output_screen("Calculating map alm..." << std::endl);
    map2alm(maskedMap_, alm, weight);
    output_screen("OK" << std::endl);

    c_.clear();
    for(int l = 0; l <= lMax_; ++l)
        c_.push_back(ps(alm, l));

    calculatePS();
}

void
Master::calculate(const Healpix_Map<double>& map)
{
    map_ = &map;

    calculate();
}

void Master::calculate(const char* mapName)
{
    Healpix_Map<double> map;
    output_screen("Reading the map..." << std::endl);
    read_Healpix_map_from_fits(std::string(mapName), map);

    if(map.Scheme() == NEST)
        map.swap_scheme();
    output_screen("OK" << std::endl);

    check(map.Scheme() == RING, "");
    
    calculate(map);
}


class Wigner3jTreats
{
public:
    Wigner3jTreats(const std::vector<double>& w, int lMax) : w_(w), coupling_(lMax), lMax_(lMax)
    {
        check(lMax_ >= 0, "");
        coupling_.resize(lMax_ + 1);
        for(int i = 0; i <= lMax_; ++i)
            coupling_[i].resize(lMax_ + 1, 0);
    }

    void process(int l1, int l2, int l3, double val)
    {
        if(l1 > lMax_ || l2 > lMax_ || l3 > lMax_)
            return;

        check(l1 >= 0, "invalid l1 = " << l1);
        check(l2 >= 0, "invalid l2 = " << l2);
        check(l3 >= 0, "invalid l3 = " << l3);

        coupling_[l1][l2] += (double(2 * l2 + 1) / (4 * Math::pi)) * double(2 * l3 + 1) * w_[l3] * val * val;
    }

    double getCoupling(int l1, int l2) const
    {
        check(l1 >= 0 && l1 <= lMax_, "invalid l1 = " << l1);
        check(l2 >= 0 && l2 <= lMax_, "invalid l2 = " << l2);
        return coupling_[l1][l2];
    }

private:
    const std::vector<double>& w_;
    std::vector<std::vector<double> > coupling_;
    const int lMax_;
};

void
Master::calculateCouplingKernel(const std::vector<double>& w, int lMax, const char* fileName)
{
    output_screen("Calculating coupling kernel..." << std::endl);
    check(w.size() >= lMax + 1, "");
    check(lMax > 0, "");

    Wigner3jTreats t(w, lMax);
    Math::Wigner3JZeroM<Wigner3jTreats> w3j(t);
    w3j.calculate(lMax * 3 / 2);
    output_screen("OK" << std::endl);

    output_screen("Saving into file " << fileName << "..." << std::endl);
    std::ofstream out(fileName);
    if(!out)
    {
        StandardException exc;
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    out << lMax << std::endl;
    for(int l1 = 0; l1 <= lMax; ++l1)
    {
        out << t.getCoupling(l1, 0);
        for(int l2 = 1; l2 <= lMax; ++l2)
            out << '\t' << t.getCoupling(l1, l2);
        out << std::endl;
    }
    out.close();
    output_screen("OK" << std::endl);
}

void
Master::calculateCouplingKernel(const Healpix_Map<double>& mask, int lMax, const char* fileName)
{
    Alm<xcomplex<double> > almMask(lMax, lMax);
    arr<double> weight(2 * mask.Nside(), 1);
    
    output_screen("Calculating alm..." << std::endl);
    map2alm(mask, almMask, weight);
    output_screen("OK" << std::endl);

    std::vector<double> w;

    for(int l = 0; l <= lMax; ++l)
        w.push_back(ps(almMask, l));

    calculateCouplingKernel(w, lMax, fileName);
}

void
Master::calculateCouplingKernel(const char* maskName, int lMax, const char* fileName)
{
    Healpix_Map<double> mask;
    output_screen("Reading the mask..." << std::endl);
    read_Healpix_map_from_fits(std::string(maskName), mask);
    
    if(mask.Scheme() == NEST)
        mask.swap_scheme();
    
    check(mask.Scheme() == RING, "");
    
    output_screen("OK" << std::endl);
    calculateCouplingKernel(mask, lMax, fileName);
}

void
Master::calculateCoupling()
{
    if(couplingKernelFileName_ != "NA")
    {
        output_screen("Reading the coupling kernel from file " << couplingKernelFileName_ << "..." << std::endl);
        std::ifstream in(couplingKernelFileName_.c_str());
        StandardException exc;
        if(!in)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Cannot open input file " << couplingKernelFileName_ << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }
        int lMax;
        in >> lMax;
        if(lMax != lMax_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Incorrect lMax = " << lMax << " in file " << couplingKernelFileName_ << ", need lMax = " << lMax_ << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }

        coupling_.resize(lMax_ + 1);
        for(int l1 = 0; l1 <= lMax_; ++l1)
        {
            coupling_[l1].resize(lMax_ + 1);
            for(int l2 = 0; l2 <= lMax_; ++l2)
                in >> coupling_[l1][l2];
        }
        in.close();
        output_screen("OK" << std::endl);
        return;
    }

    output_screen("Calculating coupling kernel..." << std::endl);
    check(coupling_.empty(), "");
    check(w_.size() == lMax_ + 1, "");

    Wigner3jTreats t(w_, lMax_);
    Math::Wigner3JZeroM<Wigner3jTreats> w3j(t);
    w3j.calculate(lMax_ * 3 / 2);
    
    coupling_.resize(lMax_ + 1);
    
    for(int l1 = 0; l1 <= lMax_; ++l1)
    {
        coupling_[l1].resize(lMax_ + 1, 0);
        for(int l2 = 0; l2 <= lMax_; ++l2)
            coupling_[l1][l2] = t.getCoupling(l1, l2);
    }
    
    /*
    std::ofstream out("coupling_kernel.txt");
    for(int l1 = 0; l1 <= lMax_; ++l1)
    {
        out << coupling_[l1][0];
        for(int l2 = 1; l2 <= lMax_; ++l2)
            out << ' ' << coupling_[l1][l2];
        out << std::endl;
    }
    out.close();
    */
    output_screen("OK" << std::endl);
}

void
Master::calculateK()
{
    output_screen("Calculating K..." << std::endl);
    int size;
    if(!bins_.empty())
        size = bins_.size() - 1;
    else
    {
        check(lMax_ > 0, "invalid lMax");
        size = lMax_ + 1;
    }
    k_.resize(size, size);
    
    if(!bins_.empty())
    {
        ProgressMeter meter(bins_.size() - 1);
        for(int b = 0; b < size; ++b)
        {
#pragma omp parallel for default(shared)
            for(int b1 = 0; b1 < size; ++b1)
            {
                k_(b, b1) = 0;
                for(int l = 0; l <= lMax_; ++l)
                {
                    for(int l1 = 0; l1 <= lMax_; ++l1)
                    {
                        k_(b, b1) += p(b, l) * coupling_[l][l1] * beam_[l1] * beam_[l1] * q(l1, b1);
                    }
                }
            }
            meter.advance();
        }
    }
    else
    {
        ProgressMeter meter(lMax_ + 1);
        for(int l = 0; l <= lMax_; ++l)
        {
            for(int l1 = 0; l1 <= lMax_; ++l1)
            {
                k_(l, l1) = coupling_[l][l1] * beam_[l1] * beam_[l1] * lFactor(l) / lFactor(l1);
            }
            meter.advance();
        }
    }
    output_screen("OK" << std::endl);
    
    //std::ofstream out("k_matrix.txt");
    //out << k_ << std::endl << std::endl;
    
    output_screen("Taking the inverse of K..." << std::endl);
    kInv_.copy(k_);
    kInv_.invert();
    //out << kInv_ << std::endl;
    //out.close();
    output_screen("OK" << std::endl);
}

void
Master::calculatePS()
{
    output_screen("Calculating final power spectrum..." << std::endl);
    int size;
    if(!bins_.empty())
        size = bins_.size() - 1;
    else
    {
        check(lMax_ > 0, "invalid lMax");
        size = lMax_ + 1;
    }
    
    ProgressMeter meter(size);
    for(int b = 0; b < size; ++b)
    {
        const double avg = (bins_.empty() ? double(b) : (bins_[b] + bins_[b + 1] - 1) / 2);
        double c = 0;
        for(int b1 = 0; b1 < size; ++b1)
        {
            if(!bins_.empty())
            {
                for(int l = 0; l <= lMax_; ++l)
                    c += kInv_(b, b1) * p(b1, l) * c_[l];
            }
            else
            {
                c += kInv_(b, b1) * c_[b1] * lFactor(b1);
            }
        }
        ps_[avg] = c;
        meter.advance();
    }
    output_screen("OK" << std::endl);
}

double
Master::p(int b, int l) const
{
    check(bins_.size() > 1, "");
    check(l >= 0 && l <= lMax_, "invalid l");
    check(b >= 0 && b < bins_.size() - 1, "invalid b");
    
    if(l < bins_[b] || l >= bins_[b + 1])
        return 0;
    
    return lFactor(l) / (bins_[b + 1] - bins_[b]);
}

double
Master::q(int l, int b) const
{
    check(bins_.size() > 1, "");
    check(l >= 0 && l <= lMax_, "invalid l");
    check(b >= 0 && b < bins_.size() - 1, "invalid b");
    
    if(l < bins_[b] || l >= bins_[b + 1])
        return 0;
    
    return 1.0 / lFactor(l);
}
