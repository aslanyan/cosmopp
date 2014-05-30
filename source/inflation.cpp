#include <cmath>
#include <limits>
#include <fstream>
#include <vector>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <math_constants.hpp>
#include <progress_meter.hpp>
#include <unit_conversions.hpp>
#include <inflation.hpp>

#include <omp.h>

SingleFieldInflation::SingleFieldInflation(const Math::RealFunction& potential, const Math::RealFunction& potentialDeriv, double phi0, double phiDot0, double phiEnd, double kPivot, double NPivot, double a0, int curvature) : K_(curvature)
{
    check(a0 > 0, "invalid starting value of a = " << a0 << ", needs to be positive");
    check(curvature == 0 || curvature == -1 || curvature == 1, "invalid curvature = " << curvature << ", needs to be 0 or +-1");
    check(phiEnd != phi0, "phi_end and phi_0 should be different");

    check(kPivot > 0, "invalid pivot point k_pivot = " << kPivot << ", must be positive");
    
    if(NPivot == 0)
    {
        NPivot = 55; // to be done better
    }

    check(NPivot > 0, "invalid N_pivot = " << NPivot << ", must be positive");

    output_screen("Evaluating the background evolution..." << std::endl);
    double t = 0;
    double a = a0;
    double phi = phi0;
    double phiDot = phiDot0;
    double logA = std::log(a);

    double zPrev = std::numeric_limits<double>::min();
    double zPrimePrev = std::numeric_limits<double>::min();

    double aPrimePrev = std::numeric_limits<double>::min();

    const double deltaLogA = 1e-3;
    
    while((phi - phiEnd) * (phi0 - phiEnd) > 0)
    {
        const double V = potential.evaluate(phi);
        const double VPrime = potentialDeriv.evaluate(phi);
        if(std::abs(V) > 1)
        {
            output_screen("WARNING: The potential energy is beyond the planck scale!!!" << std::endl << "   V = " << V << " at phi = " << phi << std::endl);
        }
        const double Hsq = - K_ / (a * a) + 1.0 / 3.0 * (phiDot * phiDot / 2.0 + V);
        check(Hsq >= 0, "total energy turned negative! H squared = " << Hsq);
        const double H = std::sqrt(Hsq);

        if(H == 0)
        {
            output_screen("Ending inflation since H = 0");
            break;
        }

        const double deltaT = deltaLogA / H;

        const double z = a * phiDot / H;
        const double zPrime = (zPrev == std::numeric_limits<double>::min() ? 0 : a * (z - zPrev) / deltaT);
        const double zDoublePrime = (zPrimePrev == std::numeric_limits<double>::min() ? 0 : a * (zPrime - zPrimePrev) / deltaT);

        const double aPrime = H * a * a;
        const double aDoublePrime = (aPrimePrev == std::numeric_limits<double>::min() ? 0 : a * (aPrime - aPrimePrev) / deltaT);

        phi_[t] = phi;
        phiDot_[t] = phiDot;
        a_[t] = a;
        H_[t] = H;
        logA_[t] = logA;
        logAInv_[logA] = t;
        aHInverse_[a * H] = t;
        z_[t] = z;
        zPrime_[t] = zPrime;
        zDoublePrime_[t] = zDoublePrime;

        aPrime_[t] = aPrime;
        aDoublePrime_[t] = aDoublePrime;

        const double deltaPhiDot = (-3 * H * phiDot - VPrime) * deltaT;
        const double deltaPhi = phiDot * deltaT;

        t += deltaT;
        phiDot += deltaPhiDot;
        phi += deltaPhi;
        logA += deltaLogA;
        a = std::exp(logA);
        
        if(zPrev != std::numeric_limits<double>::min())
            zPrimePrev = zPrime;
        zPrev = z;

        aPrimePrev = aPrime;

        if(a_.size() >= 1000000000)
        {
            StandardException exc;
            std::string exceptionStr = "Single field inflation background evolution is not converging. Stopping after 1000000000 points.";
            exc.set(exceptionStr);
            throw exc;
        }
    }

    output_screen("OK" << std::endl << "Background evaluated at " << a_.size() << " points." << std::endl);
    output_screen("Evaluating tau..." << std::endl);
    check(a_.size() > 2, "not enough points");
    Math::TableFunction<double, double>::const_iterator it = a_.end();
    --it;
    t = (*it).first;
    double tau = 0;
    tau_[t] = tau;
    while(true)
    {
        --it;
        const double tNew = (*it).first;
        a = (*it).second;
        check(a > 0, "");
        const double deltaTau = (tNew - t) / a;
        tau += deltaTau;
        t = tNew;
        tau_[t] = tau;
        if(it == a_.begin())
            break;
    }

    tauStart_ = tau;

    check(tau_.size() == a_.size(), "");

    for(it = tau_.begin(); it != tau_.end(); ++it)
    {
        t = (*it).first;
        tau = (*it).second;
        zRatioTau_[tau] = zDoublePrime_.evaluate(t) / z_.evaluate(t);
        aRatioTau_[tau] = aDoublePrime_.evaluate(t) / a_.evaluate(t);
    }
    check(zRatioTau_.size() == a_.size(), "");
    check(aRatioTau_.size() == a_.size(), "");
    output_screen("OK" << std::endl);

    it = aHInverse_.begin();
    aHStart_ = (*it).first;
    tStart_ = (*it).second;
    it = aHInverse_.end();
    --it;
    aHEnd_ = (*it).first;

    it = logA_.end();
    --it;
    const double logAEnd = (*it).second;
    check(logAEnd > NPivot, "Only " << logAEnd << " efolds during inflation, N_pivot = " << NPivot << " is too large");
    const double tPivot = logAInv_.evaluate(logAEnd - NPivot);

    phiPivot_ = phi_.evaluate(tPivot);

    it = a_.end();
    --it;
    aEnd_ = (*it).second;
    it = H_.end();
    --it;
    HEnd_ = (*it).second;

    kPivot = kUnitlessFromMpcInv(kPivot);
    aNow_ = a_.evaluate(tPivot) * H_.evaluate(tPivot) / kPivot;
}

double
SingleFieldInflation::kUnitlessFromMpcInv(double k) const
{
    check(k > 0, "");
    return 1.0 / Phys::mToUnitless(Phys::MpcToM(1.0 / k));
}

double
SingleFieldInflation::scalarP(double k) const
{
    k = kUnitlessFromMpcInv(k) * aNow_;
    check(k > 10 * aHStart_, "the mode k = " << k << " is not deep enough inside the horizon when inflation starts");
    check(k < 101 * aHEnd_, "the mode k = " << k << " is not far enough outside the horizon when inflation ends");

    double aH0 = k / 40;
    double t0 = tStart_;
    if(aH0 > aHStart_)
        t0 = aHInverse_.evaluate(aH0);

    const double tEnd = aHInverse_.evaluate(100 * k);

    const double tau0 = tau_.evaluate(t0), tauEnd = tau_.evaluate(tEnd);

    double tau = tau0;
    double v1 = std::cos(k * tau) / std::sqrt(2 * k);
    double v1Prime = -std::sqrt(k / 2) * std::sin(k * tau);
    double v2 = -std::sin(k * tau) / std::sqrt(2 * k);
    double v2Prime = -std::sqrt(k / 2) * std::cos(k * tau);

    const double deltaTauFrac = 1e-4;

    while(tau < tauEnd)
    {
        const double omega2 = k * k - zRatioTau_.evaluate(tau);
        const double deltaTau = std::min(1e-5 / k, std::abs(deltaTauFrac * tau));
        
        const double deltaV1 = v1Prime * deltaTau;
        const double deltaV2 = v2Prime * deltaTau;
        const double deltaV1Prime = -omega2 * v1 * deltaTau;
        const double deltaV2Prime = -omega2 * v2 * deltaTau;

        tau += deltaTau;
        v1 += deltaV1;
        v2 += deltaV2;
        v1Prime += deltaV1Prime;
        v2Prime += deltaV2Prime;
    }

    const double absV = v1 * v1 + v2 * v2;
    const double zEnd = z_.evaluate(tEnd);
    
    return k * k * k / (2 * Math::pi * Math::pi) * absV / (zEnd * zEnd);
}

Math::TableFunction<double, double>*
SingleFieldInflation::scalarPs(double kMin, double kMax, unsigned int nPoints) const
{
    check(kMax > kMin, "");
    check(aNow_ * kUnitlessFromMpcInv(kMax) < 101 * aHEnd_, "the maximum mode k_max = " << aNow_ * kUnitlessFromMpcInv(kMax) << " is not far enough outside the horizon when inflation ends, needs to be < " << 101 * aHEnd_);
    check(aNow_ * kUnitlessFromMpcInv(kMin) > 10 * aHStart_, "the miminum mode k_min = " << aNow_ * kUnitlessFromMpcInv(kMin) << " is not deep enough inside the horizon when inflation starts, needs to be > " << 10 * aHStart_);
    check(nPoints > 2, "");

    const double logKMin = std::log(kMin), logKMax = std::log(kMax);
    const double deltaLogK = (logKMax - logKMin) / (nPoints - 1);

    std::vector<double> res(nPoints);

    output_screen("Calculating the scalar power spectrum at " << nPoints << " points..." << std::endl);
    ProgressMeter meter(nPoints);
    omp_lock_t lock;
#pragma omp parallel for default(shared)
    for(unsigned int i = 0; i < nPoints; ++i)
    {
        double k = std::exp(logKMin + i * deltaLogK);
        if(i == 0)
            k = kMin;

        if(i == nPoints - 1)
            k = kMax;

        res[i] = scalarP(k);

        omp_set_lock(&lock);
        meter.advance();
        omp_unset_lock(&lock);
    }
    output_screen("OK" << std::endl);

    Math::TableFunction<double, double>* r = new Math::TableFunction<double, double>;
    for(unsigned int i = 0; i < nPoints; ++i)
    {
        double k = std::exp(logKMin + i * deltaLogK);
        if(i == 0)
            k = kMin;

        if(i == nPoints - 1)
            k = kMax;

        (*r)[k] = res[i];
    }

    return r;
}

double
SingleFieldInflation::tensorP(double k) const
{
    k = kUnitlessFromMpcInv(k) * aNow_;
    check(k > 10 * aHStart_, "the mode k = " << k << " is not deep enough inside the horizon when inflation starts");
    check(k < 101 * aHEnd_, "the mode k = " << k << " is not far enough outside the horizon when inflation ends");

    double aH0 = k / 40;
    double t0 = tStart_;
    if(aH0 > aHStart_)
        t0 = aHInverse_.evaluate(aH0);

    const double tEnd = aHInverse_.evaluate(100 * k);

    const double tau0 = tau_.evaluate(t0), tauEnd = tau_.evaluate(tEnd);

    double tau = tau0;
    double v1 = std::cos(k * tau) / std::sqrt(2 * k);
    double v1Prime = -std::sqrt(k / 2) * std::sin(k * tau);
    double v2 = -std::sin(k * tau) / std::sqrt(2 * k);
    double v2Prime = -std::sqrt(k / 2) * std::cos(k * tau);

    const double deltaTauFrac = 1e-4;

    while(tau < tauEnd)
    {
        const double omega2 = k * k - aRatioTau_.evaluate(tau);
        const double deltaTau = std::min(1e-5 / k, std::abs(deltaTauFrac * tau));
        
        const double deltaV1 = v1Prime * deltaTau;
        const double deltaV2 = v2Prime * deltaTau;
        const double deltaV1Prime = -omega2 * v1 * deltaTau;
        const double deltaV2Prime = -omega2 * v2 * deltaTau;

        tau += deltaTau;
        v1 += deltaV1;
        v2 += deltaV2;
        v1Prime += deltaV1Prime;
        v2Prime += deltaV2Prime;
    }

    const double absV = v1 * v1 + v2 * v2;
    const double aEnd = a_.evaluate(tEnd);
    
    return 4 * k * k * k / (Math::pi * Math::pi) * absV / (aEnd * aEnd);
}

Math::TableFunction<double, double>*
SingleFieldInflation::tensorPs(double kMin, double kMax, unsigned int nPoints) const
{
    check(kMax > kMin, "");
    check(aNow_ * kUnitlessFromMpcInv(kMax) < 101 * aHEnd_, "the maximum mode k_max = " << aNow_ * kUnitlessFromMpcInv(kMax) << " is not far enough outside the horizon when inflation ends, needs to be < " << 101 * aHEnd_);
    check(aNow_ * kUnitlessFromMpcInv(kMin) > 10 * aHStart_, "the miminum mode k_min = " << aNow_ * kUnitlessFromMpcInv(kMin) << " is not deep enough inside the horizon when inflation starts, needs to be > " << 10 * aHStart_);
    check(nPoints > 2, "");

    const double logKMin = std::log(kMin), logKMax = std::log(kMax);
    const double deltaLogK = (logKMax - logKMin) / (nPoints - 1);

    std::vector<double> res(nPoints);

    output_screen("Calculating the tensor power spectrum at " << nPoints << " points..." << std::endl);
    ProgressMeter meter(nPoints);
    omp_lock_t lock;
#pragma omp parallel for default(shared)
    for(unsigned int i = 0; i < nPoints; ++i)
    {
        double k = std::exp(logKMin + i * deltaLogK);
        if(i == 0)
            k = kMin;

        if(i == nPoints - 1)
            k = kMax;

        res[i] = tensorP(k);

        omp_set_lock(&lock);
        meter.advance();
        omp_unset_lock(&lock);
    }
    output_screen("OK" << std::endl);

    Math::TableFunction<double, double>* r = new Math::TableFunction<double, double>;
    for(unsigned int i = 0; i < nPoints; ++i)
    {
        double k = std::exp(logKMin + i * deltaLogK);
        if(i == 0)
            k = kMin;

        if(i == nPoints - 1)
            k = kMax;

        (*r)[k] = res[i];
    }

    return r;
}
