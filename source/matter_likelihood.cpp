#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <scale_factor.hpp>
#include <unit_conversions.hpp>
#include <matter_likelihood.hpp>

#include <gmd.h>
#include <lavd.h>
#include <laslv.h>
#include <lavli.h>
#include <blas2pp.h>
#include <blas3pp.h>

MatterLikelihood::MatterLikelihood(const char* pkFileName, const char* covFileName, bool isLog, double khMin, double khMax) : scale_(false)
{
    readPk(pkFileName, khMin, khMax);
    if(isLog)
        readLogCov(covFileName);
    else
        readCov(covFileName);
}

double
MatterLikelihood::calculateLin(const Math::RealFunction& matterPk, const CosmologicalParams& params) const
{
    check(!data_.empty(), "need to initialize the data first");

    const int n = data_.size();

    check(kh_.size() == n, "");
    check(cInv_.size(0) == n && cInv_.size(1) == n, "need to initialize the covariance matrix first");

    const double h = params.getH();

    double alpha = 1;
    if(scale_)
    {
        const double dv = DV(params, z_);
        alpha = dvFid_ / dv;
    }

    std::vector<double> p1(n), p2(n);
    for(int i = 0; i < n; ++i)
    {
        const double kh = kh_[i] * alpha;
        const double k = kh * h;
        const double p = matterPk.evaluate(k) / (alpha * alpha * alpha);

        p1[i] = p * h * h * h;
        p2[i] = 1.0;
    }

    LaGenMatDouble A(2, 2);

    A(0, 0) = 0;
    A(0, 1) = 0;
    A(1, 0) = 0;
    A(1, 1) = 0;

    double dataLike = 0;
    double b0 = 0, b1 = 0;

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            const double c = cInv_(i, j);
            const double p1c = p1[i] * c, p2c = p2[i] * c;
            const double dj = data_[j];

            A(0, 0) += p1c * p1[j];
            A(1, 0) += p1c * p2[j];
            A(1, 1) += p2c * p2[j];

            b0 += p1c * dj;
            b1 += p2c * dj;

            dataLike += data_[i] * c * dj;
        }
    }

    A(0, 1) = A(1, 0);

    const double detA = A(0, 0) * A(1, 1) - A(1, 0) * A(1, 0);

    // inverting A
    LaVectorLongInt pivotA(2);
    LUFactorizeIP(A, pivotA);
    LaLUInverseIP(A, pivotA);

    const double bAb = b0 * A(0, 0) * b0 + b0 * A(0, 1) * b1 + b1 * A(1, 0) * b0 + b1 * A(1, 1) * b1;

    const double like = std::log(detA) - bAb + dataLike;

    return like;
    /*
    check(!data_.empty(), "need to initialize the data first");
    
    const int n = data_.size();
    check(kh_.size() == n, "");
    check(cInv_.size(0) == n && cInv_.size(1) == n, "need to initialize the covariance matrix first");

    const double h = params.getH();

    double alpha = 1;
    if(scale_)
    {
        const double dv = DV(params, z_);
        alpha = dvFid_ / dv;
    }

    std::vector<double> p1(n);
    for(int i = 0; i < n; ++i)
    {
        const double kh = kh_[i] * alpha;
        const double k = kh * h;
        const double p = matterPk.evaluate(k) / (alpha * alpha * alpha);

        p1[i] = p * h * h * h;
    }

    // marginalize over b2
    if(b2 == 0)
    {
        double dCInvd = 0, dCInvp = 0, pCInvp = 0;

        for(int i = 0; i < n; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                dCInvd += data_[i] * cInv_(i, j) * data_[j];
                dCInvp += data_[i] * cInv_(i, j) * p1[j];
                pCInvp += p1[i] * cInv_(i, j) * p1[j];
            }
        }

        return dCInvd - dCInvp * dCInvp / pCInvp + std::log(pCInvp);
    }

    double res = 0;
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
            res += (b2 * p1[i] - data_[i]) * cInv_(i, j) * (b2 * p1[j] - data_[j]);
    }

    return res;
    */
}

double
MatterLikelihood::calculate(const Math::RealFunction& matterPk, const CosmologicalParams& params) const
{
    check(!data_.empty(), "need to initialize the data first");

    const int n = data_.size();

    check(kh_.size() == n, "");
    check(cInv_.size(0) == n && cInv_.size(1) == n, "need to initialize the covariance matrix first");

    const double h = params.getH();

    double alpha = 1;
    if(scale_)
    {
        const double dv = DV(params, z_);
        alpha = dvFid_ / dv;
    }

    std::vector<double> p1(n), p2(n);
    for(int i = 0; i < n; ++i)
    {
        const double kh = kh_[i] * alpha;
        const double k = kh * h;
        const double p = matterPk.evaluate(k) / (alpha * alpha * alpha);

        p1[i] = p * h * h * h / (1.0 + 1.4 * kh);
        p2[i] = p1[i] * kh * kh;
    }

    LaGenMatDouble A(2, 2);

    A(0, 0) = 0;
    A(0, 1) = 0;
    A(1, 0) = 0;
    A(1, 1) = 0;

    double dataLike = 0;
    double b0 = 0, b1 = 0;

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            const double c = cInv_(i, j);
            const double p1c = p1[i] * c, p2c = p2[i] * c;
            const double dj = data_[j];

            A(0, 0) += p1c * p1[j];
            A(1, 0) += p1c * p2[j];
            A(1, 1) += p2c * p2[j];

            b0 += p1c * dj;
            b1 += p2c * dj;

            dataLike += data_[i] * c * dj;
        }
    }

    A(0, 1) = A(1, 0);

    const double detA = A(0, 0) * A(1, 1) - A(1, 0) * A(1, 0);

    // inverting A
    LaVectorLongInt pivotA(2);
    LUFactorizeIP(A, pivotA);
    LaLUInverseIP(A, pivotA);

    const double bAb = b0 * A(0, 0) * b0 + b0 * A(0, 1) * b1 + b1 * A(1, 0) * b0 + b1 * A(1, 1) * b1;

    const double like = std::log(detA) - bAb + dataLike;

    return like;
}

void
MatterLikelihood::useScaling(const CosmologicalParams& paramsFid, double z)
{
    scale_ = true;
    z_ = z;
    dvFid_ = DV(paramsFid, z);
}

double
MatterLikelihood::DV(const CosmologicalParams& params, double z) const
{
    check(z >= 0, "z cannot be negative");

    ScaleFactorFunctionClass sf;
    sf.initialize(params);

    const double age = sf.age();
    const double t1 = sf.time(1.0 / (1.0 + z));

    const double com = sf.comovingDistance(t1, age);

    const double comMpc = Phys::mToMpc(Phys::unitlessToM(com));
    const double hz = sf.hubble(1.0 / (1.0 + z));
    const double invHMpc = Phys::mToMpc(Phys::unitlessToM(1.0 / hz));

    const double dv = std::pow(z * comMpc * comMpc * invHMpc, 1.0 / 3.0);
    return dv;
}

void
MatterLikelihood::readPk(const char* fileName, double khMin, double khMax)
{
    StandardException exc;
    std::ifstream in(fileName);

    if(!in)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot open input file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    kh_.clear();
    data_.clear();

    nIgnored_ = 0;
    nTotal_ = 0;

    while(!in.eof())
    {
        std::string s;
        std::getline(in, s);

        if(s[0] == '#')
            continue;

        if(s == "")
            break;

        std::stringstream str(s);
        double kh, pk;

        str >> kh >> pk;

        ++nTotal_;

        if(kh < khMin)
        {
            ++nIgnored_;
            continue;
        }
        if(kh > khMax)
            continue;

        kh_.push_back(kh);
        data_.push_back(pk);
    }

    in.close();
}

void
MatterLikelihood::readCov(const char* fileName)
{
    check(!data_.empty(), "need to initialize the data first");
    
    const int n = data_.size();
    check(kh_.size() == n, "");

    StandardException exc;
    std::ifstream in(fileName);

    if(!in)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot open input file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    cInv_.resize(n, n);

    std::string s;
    for(int i = 0; i < nIgnored_; ++i)
        std::getline(in, s);

    for(int i = 0; i < n; ++i)
    {
        std::getline(in, s);
        double x;
        std::stringstream str(s);
        for(int j = 0; j < nIgnored_; ++j)
            str >> x;

        for(int j = 0; j < n; ++j)
        {
            str >> x;
            cInv_(i, j) = x;
        }
    }
    
    in.close();

    std::ofstream out("matter_test.txt");
    for(int i = 0; i < n; ++i)
        out << kh_[i] << '\t' << data_[i] << '\t' << std::sqrt(cInv_(i, i)) << std::endl;
    out.close();

    // inverting the covariance matrix
    LaVectorLongInt pivotC(n);
    LUFactorizeIP(cInv_, pivotC);
    LaLUInverseIP(cInv_, pivotC);
}

void
MatterLikelihood::readLogCov(const char* fileName)
{
    check(!data_.empty(), "need to initialize the data first");
    
    const int n = data_.size();
    check(kh_.size() == n, "");

    StandardException exc;
    std::ifstream in(fileName);

    if(!in)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot open input file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    cInv_.resize(n, n);

    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            cInv_(i, j) = 0;

    while(!in.eof())
    {
        std::string s;
        std::getline(in, s);

        if(s[0] == '#')
            continue;

        if(s == "")
            break;

        std::stringstream str(s);

        int i, j;
        double c;
        str >> i >> j >> c;

        if(i < 1 || i > nTotal_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Invalid format of the input file " << fileName << ": the first index is " << i << " but needs to be between 1 and " << nTotal_ << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }

        if(j < 1 || j > nTotal_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Invalid format of the input file " << fileName << ": the second index is " << j << " but needs to be between 1 and " << nTotal_ << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }

        c *= (std::log(10) * std::log(10));

        --i;
        --j;

        i -= nIgnored_;
        j -= nIgnored_;
        
        if(i >= 0 && i < n && j >= 0 && j < n)
            cInv_(i, j) = data_[i] * data_[j] * (c > 0.001 ? (std::exp(c) - 1.0) : c);
    }
    
    in.close();

    std::ofstream out("matter_test.txt");
    for(int i = 0; i < n; ++i)
        out << kh_[i] << '\t' << data_[i] << '\t' << std::sqrt(cInv_(i, i)) << std::endl;
    out.close();

    // inverting the covariance matrix
    LaVectorLongInt pivotC(n);
    LUFactorizeIP(cInv_, pivotC);
    LaLUInverseIP(cInv_, pivotC);
}
