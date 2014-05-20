#include <cmath>
#include <limits>
#include <fstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <math_constants.hpp>
#include <inflation.hpp>

SingleFieldInflation::SingleFieldInflation(const Math::RealFunction& potential, const Math::RealFunction& potentialDeriv, double phi0, double phiDot0, double phiEnd, double a0, int curvature) : K_(curvature)
{
    check(a0 > 0, "invalid starting value of a = " << a0 << ", needs to be positive");
    check(curvature == 0 || curvature == -1 || curvature == 1, "invalid curvature = " << curvature << ", needs to be 0 or +-1");
    check(phiEnd != phi0, "phi_end and phi_0 should be different");

    output_screen("Evaluating the background evolution..." << std::endl);
    double t = 0;
    double a = a0;
    double phi = phi0;
    double phiDot = phiDot0;
    double logA = std::log(a);

    double zPrev = std::numeric_limits<double>::min();
    double zPrimePrev = std::numeric_limits<double>::min();

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

        phi_[t] = phi;
        phiDot_[t] = phiDot;
        a_[t] = a;
        H_[t] = H;
        logA_[t] = logA;
        aHInverse_[a * H] = t;
        z_[t] = z;
        zPrime_[t] = zPrime;
        zDoublePrime_[t] = zDoublePrime;

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
    }
    check(zRatioTau_.size() == a_.size(), "");

    it = aHInverse_.begin();
    aHStart_ = (*it).first;
    tStart_ = (*it).second;
    it = aHInverse_.end();
    --it;
    aHEnd_ = (*it).first;
    output_screen("OK" << std::endl);
}

double
SingleFieldInflation::scalarP(double k) const
{
    check(k > 101 * aHStart_, "the mode k = " << k << " is not deep enough inside the horizon when inflation starts");
    check(k < 1001 * aHEnd_, "the mode k = " << k << " is not far enough outside the horizon when inflation ends");

    double aH0 = k / 1000;
    double t0 = tStart_;
    if(aH0 > aHStart_)
        t0 = aHInverse_.evaluate(aH0);

    const double tEnd = aHInverse_.evaluate(1000 * k);

    const double tau0 = tau_.evaluate(t0), tauEnd = tau_.evaluate(tEnd);

    //Math::TableFunction<double, double> vRe, vPrimeRe, vIm, vPrimeIm;
    
    double tau = tau0;
    double v1 = std::cos(k * tau) / std::sqrt(2 * k);
    double v1Prime = -std::sqrt(k / 2) * std::sin(k * tau);
    double v2 = -std::sin(k * tau) / std::sqrt(2 * k);
    double v2Prime = -std::sqrt(k / 2) * std::cos(k * tau);

    const double deltaTauFrac = 1e-3;

    std::ofstream outTest("test_mode.txt");
    while(tau < tauEnd)
    {
        //vRe[tau] = v1;
        //vPrimeRe[tau] = v1Prime;
        //vIm[tau] = v2;
        //vPrimeIm[tau] = v2Prime;

        outTest << tau << '\t' << v1 << '\t' << v2 << '\t' << zRatioTau_.evaluate(tau) << std::endl;

        const double omega2 = k * k - zRatioTau_.evaluate(tau);
        const double deltaTau = std::min(1e-4 / k, std::abs(deltaTauFrac * tau));
        
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

