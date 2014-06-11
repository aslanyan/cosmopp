#include <string>
#include <sstream>
#include <ctime>
#include <fstream>
#include <algorithm>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <math_constants.hpp>
#include <progress_meter.hpp>
#include <utils.hpp>
#include <conjugate_gradient.hpp>
#include <numerics.hpp>
#include <cmb_gibbs.hpp>

#include <chealpix.h>
#include <healpix_base.h>
#include <alm.h>
#include <alm_healpix_tools.h>
#include <xcomplex.h>
#include <fitshandle.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>

CMBGibbsSampler::CMBGibbsSampler(const char* mapName, const char* noiseMapName, const char* maskName, double pixelNoise, int lMax, double fwhm, const char* startingClFileName, time_t seed)
{
    Healpix_Map<double> map, noise, mask;
    read_Healpix_map_from_fits(std::string(mapName), map);
    read_Healpix_map_from_fits(std::string(noiseMapName), noise);
    read_Healpix_map_from_fits(std::string(maskName), mask);

    std::vector<double> cl;

    StandardException exc;
    std::ifstream in(startingClFileName);
    if(!in)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot open input file " << startingClFileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    while(!in.eof())
    {
        std::string s;
        std::getline(in, s);
        if(s == "")
            break;
        
        std::stringstream str(s);
        double c;
        str >> c;
        cl.push_back(c);
    }
    
    in.close();

    CMBGibbsSampler(map, noise, mask, pixelNoise, lMax, fwhm, cl, seed);
}

CMBGibbsSampler::CMBGibbsSampler(const Healpix_Map<double>& map, const Healpix_Map<double>& noiseMap, const Healpix_Map<double>& mask, double pixelNoise, int lMax, double fwhm, const std::vector<double>& startingCl, time_t seed)
{
    StandardException exc;
    if(map.Nside() != noiseMap.Nside())
    {
        std::stringstream exceptionStr;
        exceptionStr << "The map has N_side = " << map.Nside() << " while the noise map has N_side = " << noiseMap.Nside() << ". They must be the same.";
        exc.set(exceptionStr.str());
        throw exc;
    }

    if(map.Nside() != mask.Nside())
    {
        std::stringstream exceptionStr;
        exceptionStr << "The map has N_side = " << map.Nside() << " while the mask has N_side = " << mask.Nside() << ". They must be the same.";
        exc.set(exceptionStr.str());
        throw exc;
    }

    if(map.Scheme() != noiseMap.Scheme())
    {
        std::string exceptionStr = "The map and the noise map have different ordering schemes.";
        exc.set(exceptionStr);
        throw exc;
    }

    map_.SetNside(map.Nside(), map.Scheme());
    for(long i = 0; i < map.Npix(); ++i)
        map_[i] = map[i] + noiseMap[i];

    if(map_.Scheme() == NEST)
        map_.swap_scheme();

    check(map_.Scheme() == RING, "");

    mask_.SetNside(mask.Nside(), mask.Scheme());
    for(long i = 0; i < mask.Npix(); ++i)
        mask_[i] = mask[i];

    if(mask_.Scheme() == NEST)
        mask_.swap_scheme();

    check(mask_.Scheme() == RING, "");

    check(lMax >= 2, "");
    check(startingCl.size() >= lMax + 1, "");

    lMax_ = lMax;

    cl_.resize(lMax + 1);
    sigmaL_.resize(lMax + 1);

    for(int l = 0; l <= lMax; ++l)
        cl_[l] = startingCl[l];

    beam_.resize(lMax + 1, 1);
    Utils::readPixelWindowFunction(beam_, map_.Nside(), lMax, fwhm);

    check(pixelNoise > 0, "invalid pixel noise " << pixelNoise);
    pixelNoise_ = pixelNoise;

    a00_ = 0;
    a1m1_ = 0;
    a10_ = 0;
    a11_ = 0;

    y00_.SetNside(map.Nside(), RING);
    y1m1_.SetNside(map.Nside(), RING);
    y10_.SetNside(map.Nside(), RING);
    y11_.SetNside(map.Nside(), RING);

    for(long i = 0; i < y00_.Npix(); ++i)
    {
        double theta, phi;
        pix2ang_ring(y00_.Nside(), i, &theta, &phi);

        y00_[i] = 1.0 / std::sqrt(4.0 * Math::pi);
        const double factor = std::sqrt(3.0 / (4.0 * Math::pi));
        y1m1_[i] = factor * std::sin(theta) * std::sin(phi);
        y10_[i] = factor * std::cos(theta);
        y11_[i] = factor * std::sin(theta) * std::cos(phi);
    }

    Healpix_Map<double> signal;
    signal.SetNside(map.Nside(), RING);
    for(long i = 0; i < map.Npix(); ++i)
        signal[i] = map[i];

    if(signal.Scheme() == NEST)
        signal.swap_scheme();

    check(signal.Scheme() == RING, "");

    s_.Set(lMax, lMax);
    s_.SetToZero();
    arr<double> weight(2 * signal.Nside(), 1);
    map2alm(signal, s_, weight);

    for(int l = 0; l <= lMax; ++l)
        for(int m = 0; m <= l; ++m)
            s_(l, m) /= beam_[l];

    calculateSigmaL();

    if(seed == 0)
        seed = std::time(0);

    generator_ = new Math::GaussianGenerator(seed, 0, 1);
}

CMBGibbsSampler::~CMBGibbsSampler()
{
    delete generator_;
}

void
CMBGibbsSampler::generateCl()
{
    cl_[0] = 0;
    cl_[1] = 0;
    for(int l = 2; l <= lMax_; ++l)
    {
        double rho = 0;
        for(int i = 0; i < 2 * l - 1; ++i)
        {
            const double rand = generator_->generate();
            rho += rand * rand;
        }

        cl_[l] = sigmaL_[l] / rho;
    }
}

class CmbGibbsMDTreats
{
public:
    CmbGibbsMDTreats(const Healpix_Map<double>& y00, const Healpix_Map<double>& y1m1, const Healpix_Map<double>& y10, const Healpix_Map<double>& y11, const Healpix_Map<double>& mask, double pixelNoise) : y00_(y00), y1m1_(y1m1), y10_(y10), y11_(y11), mask_(mask), pixelNoise_(pixelNoise)
    {
    }

    void multiplyByMatrix(const std::vector<double>& original, std::vector<double>& result) const
    {
        check(original.size() == 4, "");
        check(result.size() == 4, "");

        Healpix_Map<double> map;
        map.SetNside(y00_.Nside(), RING);
        for(long i = 0; i < map.Npix(); ++i)
        {
            map[i] = y00_[i] * original[0] + y1m1_[i] * original[1] + y10_[i] * original[2] + y11_[i] * original[3];

            if(mask_[i] > 0.5)
                map[i] /= (pixelNoise_ * pixelNoise_);
            else
                map[i] = 0;
        }

        for(int i = 0; i < 4; ++i)
            result[i] = 0;

        for(long i = 0; i < map.Npix(); ++i)
        {
            result[0] += y00_[i] * map[i];
            result[1] += y1m1_[i] * map[i];
            result[2] += y10_[i] * map[i];
            result[3] += y11_[i] * map[i];
        }
    }

    void preconditioner(const std::vector<double>& original, std::vector<double>& result) const
    {
        result = original;
    }

private:
    const Healpix_Map<double>& y00_, y1m1_, y10_, y11_, mask_;
    double pixelNoise_;
};

void
CMBGibbsSampler::generateW()
{
    Alm<xcomplex<double> > alm(lMax_, lMax_);
    for(int l = 0; l <= lMax_; ++l)
        for(int m = 0; m <= l; ++m)
            alm(l, m) = s_(l, m) * beam_[l];

    Healpix_Map<double> map, omega0, omega1, omega2, omega3;
    map.SetNside(map_.Nside(), RING);
    omega0.SetNside(map_.Nside(), RING);
    omega1.SetNside(map_.Nside(), RING);
    omega2.SetNside(map_.Nside(), RING);
    omega3.SetNside(map_.Nside(), RING);
    alm2map(alm, map);

    std::vector<double> b(4, 0.0);

    for(long i = 0; i < map.Npix(); ++i)
    {
        map[i] = map_[i] - map[i];
        omega0[i] = generator_->generate();
        omega1[i] = generator_->generate();
        omega2[i] = generator_->generate();
        omega3[i] = generator_->generate();

        if(mask_[i] > 0.5)
        {
            map[i] /= (pixelNoise_ * pixelNoise_);
            omega0[i] /= pixelNoise_;
            omega1[i] /= pixelNoise_;
            omega2[i] /= pixelNoise_;
            omega3[i] /= pixelNoise_;
        }
        else
        {
            map[i] = 0;
            omega0[i] = 0;
            omega1[i] = 0;
            omega2[i] = 0;
            omega3[i] = 0;
        }

        map[i] += (omega0[i] + omega1[i] + omega2[i] + omega3[i]);

        b[0] += y00_[i] * map[i];
        b[1] += y1m1_[i] * map[i];
        b[2] += y10_[i] * map[i];
        b[3] += y11_[i] * map[i];
    }

    CmbGibbsMDTreats treats(y00_, y1m1_, y10_, y11_, mask_, pixelNoise_);
    Math::ConjugateGradient<CmbGibbsMDTreats> cg(4, &treats, b);
    const std::vector<double>& res = cg.solve(1e-20);
    check(res.size() == 4, "");

    a00_ = res[0];
    a1m1_ = res[1];
    a10_ = res[2];
    a11_ = res[3];
}

class CmbGibbsCGTreats
{
public:
    CmbGibbsCGTreats(long nSide, int lMax, const std::vector<double>& cl, const Healpix_Map<double>& mask, double pixelNoise, const std::vector<double>& beam) : nSide_(nSide), lMax_(lMax), cl_(cl), mask_(mask), pixelNoise_(pixelNoise), beam_(beam)
    {
        check(lMax_ >= 2, "");
        check(cl_.size() == lMax_ + 1, "");
        check(pixelNoise_ > 0, "");
        check(beam_.size() == lMax_ + 1, "");
        check(mask_.Nside() == nSide_, "");

        goodPix_ = 0;
        for(long i = 0; i < mask_.Npix(); ++i)
            if(mask_[i] > 0.5)
                ++goodPix_;
    }

    static void almToVector(const Alm<xcomplex<double> >& alm, std::vector<double>& v)
    {
        const int lMax = alm.Lmax();
        v.clear();
        for(int l = 0; l <= lMax; ++l)
        {
            v.push_back(alm(l, 0).real());
            check(Math::areEqual(alm(l, 0).imag(), 0.0, 1e-10), "");

            for(int m = 1; m <= l; ++m)
            {
                v.push_back(alm(l, m).real());
                v.push_back(alm(l, m).imag());
            }
        }
    }

    static void vectorToAlm(const std::vector<double>& v, Alm<xcomplex<double> > & alm)
    {
        const int lMax = alm.Lmax();
        std::vector<double>::const_iterator it = v.begin();
        for(int l = 0; l <= lMax; ++l)
        {
            check(it != v.end(), "");
            alm(l, 0) = xcomplex<double>((*it), 0.0);
            ++it;
            for(int m = 1; m <= l; ++m)
            {
                check(it != v.end(), "");
                const double re = *it;
                ++it;
                check(it != v.end(), "");
                const double im = *it;
                ++it;
                alm(l, m) = xcomplex<double>(re, im);
            }
        }
    }

    void multiplyByMatrix(const std::vector<double>& original, std::vector<double>& result) const
    {
        check(result.size() == original.size(), "");

        Alm<xcomplex<double> > alm(lMax_, lMax_);
        vectorToAlm(original, alm);

        Healpix_Map<double> map;
        map.SetNside(nSide_, RING);

        // perform the matrix operation, part 1
        for(int l = 0; l <= lMax_; ++l)
        {
            alm(l, 0) *= (beam_[l] * std::sqrt(cl_[l]));
            for(int m = 1; m <= l; ++m)
                alm(l, m) *= (beam_[l] * std::sqrt(cl_[l]));
        }

        // back to pixel space
        alm2map(alm, map);

        // apply noise
        for(long i = 0; i < map.Npix(); ++i)
        {
            if(mask_[i] > 0.5)
                map[i] /= (pixelNoise_ * pixelNoise_);
            else
                map[i] = 0;
        }

        // back to harmonic space
        arr<double> weight(2 * map.Nside(), map.Npix() / (4 * Math::pi));
        map2alm(map, alm, weight);

        // perform the matrix operation, part 2
        for(int l = 0; l <= lMax_; ++l)
        {
            alm(l, 0) *= (beam_[l] * std::sqrt(cl_[l]));
            for(int m = 1; m <= l; ++m)
                alm(l, m) *= (beam_[l] * std::sqrt(cl_[l]));
        }

        almToVector(alm, result);

        for(int i = 0; i < result.size(); ++i)
            result[i] += original[i];
    }

    void preconditioner(const std::vector<double>& original, std::vector<double>& result) const
    {
        // go to harmonic space
        Alm<xcomplex<double> > alm(lMax_, lMax_);
        vectorToAlm(original, alm);

        // perform preconditioning
        for(int l = 0; l <= lMax_; ++l)
            for(int m = 0; m <= l; ++m)
            {
                const double nl = pixelNoise_ * pixelNoise_ * 4 * Math::pi / goodPix_;
                const double d = 1.0 + cl_[l] * beam_[l] * beam_[l] / nl;
                alm(l, m) /= d;
            }

        // write map back into vector
        almToVector(alm, result);
    }

private:
    long nSide_;
    int lMax_;
    const std::vector<double>& cl_;
    const Healpix_Map<double>& mask_;
    double pixelNoise_;
    const std::vector<double>& beam_;
    int goodPix_;
};

void
CMBGibbsSampler::generateSignal()
{
    // start in pixel space
    Healpix_Map<double> map, omega0, omega1;
    map.SetNside(map_.Nside(), RING);
    omega0.SetNside(map_.Nside(), RING);
    omega1.SetNside(map_.Nside(), RING);

    for(long i = 0; i < map_.Npix(); ++i)
    {
        omega0[i] = generator_->generate();
        omega1[i] = generator_->generate();
        map[i] = map_[i] - a00_ * y00_[i] - a1m1_ * y1m1_[i] - a10_ * y10_[i] - a11_ * y11_[i];

        if(mask_[i] > 0.5)
        {
            map[i] /= (pixelNoise_ * pixelNoise_);
            omega1[i] /= pixelNoise_;
        }
        else
        {
            map[i] = 0;
            omega1[i] = 0;
        }
    }

    // harmonic space
    Alm<xcomplex<double> > alm(lMax_, lMax_), alm1(lMax_, lMax_), alm2(lMax_, lMax_);
    arr<double> weight(2 * map.Nside(), map.Npix() / (4 * Math::pi));
    arr<double> weight1(2 * map.Nside(), std::sqrt(map.Npix() / (4 * Math::pi)));
    map2alm(map, alm, weight);
    map2alm(omega1, alm1, weight);

    for(int l = 0; l <= lMax_; ++l)
        for(int m = 0; m <= l; ++m)
            alm(l, m) += alm1(l, m);

    // beam and sqrt(cl)
    for(int l = 0; l <= lMax_; ++l)
    {
        alm(l, 0) *= (beam_[l] * std::sqrt(cl_[l]));
        for(int m = 1; m <= l; ++m)
            alm(l, m) *= (beam_[l] * std::sqrt(cl_[l]));
    }

    arr<double> weight2(2 * map.Nside(), 1);
    map2alm(omega0, alm2, weight1);

    for(int l = 0; l <= lMax_; ++l)
        for(int m = 0; m <= l; ++m)
            alm(l, m) += alm2(l, m);
    std::vector<double> b;
    CmbGibbsCGTreats::almToVector(alm, b);

    CmbGibbsCGTreats cgTreats(map.Nside(), lMax_, cl_, mask_, pixelNoise_, beam_);

    Math::ConjugateGradient<CmbGibbsCGTreats> cg(b.size(), &cgTreats, b);

    //output_screen("Generating the signal..." << std::endl);
    int iter;
    const std::vector<double>& res = cg.solve(1e-20, &iter);
    output_screen2("The CG converged in " << iter << " iterations." << std::endl);

#ifdef CHECKS_ON
    std::vector<double> x(res.size());
    cgTreats.multiplyByMatrix(res, x);
    double diff = 0;
    for(int i = 0; i < x.size(); ++i)
        diff += (x[i] - b[i]) * (x[i] - b[i]);

    diff = std::sqrt(diff);
    check(diff < 1e-5, "");
#endif

    /*
    // into a map
    for(long i = 0; i < map.Npix(); ++i)
        map[i] = res[i];

    // harmonic space
    map2alm(map, s_, weight);
    */
    CmbGibbsCGTreats::vectorToAlm(res, s_);

    // fix for sqrt(cl)
    for(int l = 0; l <= lMax_; ++l)
    {
        s_(l, 0) *= std::sqrt(cl_[l]);
        for(int m = 1; m <= l; ++m)
            s_(l, m) *= std::sqrt(cl_[l]);
    }

    calculateSigmaL();
}

void
CMBGibbsSampler::step()
{
    generateW();
    generateSignal();
    generateCl();
}

void
CMBGibbsSampler::calculateSigmaL()
{
    for(int l = 0; l <= lMax_; ++l)
    {
        sigmaL_[l] = s_(l, 0).norm();
        for(int m = 1; m <= l; ++m)
            sigmaL_[l] += 2 * s_(l, m).norm();
    }
}

double
CMBGibbsSampler::calculateLikelihood(const std::vector<double>& cl, const std::vector<double>& sigmaL, int lMax)
{
    check(lMax >= 2, "");
    check(cl.size() >= lMax + 1, "");
    check(sigmaL.size() >= lMax + 1, "");

    double res = 0;
    for(int l = 2; l <= lMax; ++l)
    {
        check(cl[l] > 0, "");

        double sigma = sigmaL[l] / double(2 * l + 1);
        res -= (2 * l - 1) * std::log(sigma);
        res += (2 * l + 1) * std::log(cl[l]);
        res += (2 * l + 1) * sigma / cl[l];
    }

    return res;
}

void
CMBGibbsSampler::generateChain(GibbsSampleChain& chain, int nSamples, int burnIn)
{
    check(burnIn >= 0, "invalid burnIn = " << burnIn);
    check(nSamples >= 1, "invalid nSamples = " << nSamples);

    output_screen("Generating the Gibbs chain..." << std::endl);

    ProgressMeter meter(nSamples + burnIn);

    chain.clear();
    for(int i = 0; i < burnIn; ++i)
    {
        step();
        meter.advance();
    }

    for(int i = 0; i < nSamples; ++i)
    {
        step();
        chain.push_back(sigmaL_);
        meter.advance();
    }
    output_screen("OK" << std::endl);
}

double
CMBGibbsSampler::calculateLikelihood(const std::vector<double>& cl, const GibbsSampleChain& chain, int lMax)
{
    check(lMax >= 2, "");
    check(cl.size() >= lMax + 1, "");
    check(!chain.empty(), "the chain is empty");

    double like = 0;
    std::vector<double> likeVec(chain.size(), 0);

#pragma omp parallel for default(shared)
    for(int i = 0; i < chain.size(); ++i)
    {
        check(chain[i].size() >= lMax + 1, "");
        likeVec[i] = calculateLikelihood(cl, chain[i], lMax);
    }

    std::sort(likeVec.begin(), likeVec.end());

    for(int i = 0; i < likeVec.size(); ++i)
    {
        if(likeVec[i] - likeVec[0] < 20)
            like += std::exp(-(likeVec[i] - likeVec[0]) / 2.0);
    }

    like /= chain.size();
    return likeVec[0] - 2 * std::log(like);
}

void
CMBGibbsSampler::writeChainIntoFile(const GibbsSampleChain& chain, const char* fileName)
{
    check(!chain.empty(), "");
    check(!chain[0].empty(), "");

    StandardException exc;
    std::ofstream out(fileName, std::ios::binary | std::ios::out);
    if(!out)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into output file " << fileName;
        exc.set(exceptionStr.str());
        throw exc;
    }

    const int size = chain.size();
    const int lMax = chain[0].size() - 1;
    out.write((char*)(&size), sizeof(size));
    out.write((char*)(&lMax), sizeof(lMax));
    for(int i = 0; i < size; ++i)
        out.write((char*)(&(chain[i][0])), (lMax + 1) * sizeof(double));

    out.close();
}

void
CMBGibbsSampler::readChainFromFile(GibbsSampleChain& chain, const char* fileName)
{
    chain.clear();

    StandardException exc;
    std::ifstream in(fileName, std::ios::in | std::ios::binary);
    if(!in)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot read the input file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    int size, lMax;
    in.read((char*)(&size), sizeof(size));
    if(size <= 0)
    {
        std::stringstream exceptionStr;
        exceptionStr << "The input file " << fileName << " contains an invalid chain. The size of the chain is " << size << ", it needs to be positive.";
        exc.set(exceptionStr.str());
        throw exc;
    }

    in.read((char*)(&lMax), sizeof(lMax));
    if(lMax < 2)
    {
        std::stringstream exceptionStr;
        exceptionStr << "The input file " << fileName << " contains an invalid chain. l_max = " << lMax << ", it needs to be >= 2.";
        exc.set(exceptionStr.str());
        throw exc;
    }

    chain.resize(size);
    for(int i = 0; i < size; ++i)
    {
        chain[i].resize(lMax + 1);
        in.read((char*)(&(chain[i][0])), (lMax + 1) * sizeof(double));
    }

    in.close();
}
