#ifndef COSMO_PP_INFLATION_HPP
#define COSMO_PP_INFLATION_HPP

#include <function.hpp>
#include <table_function.hpp>

class SingleFieldInflation
{
public:
    SingleFieldInflation(const Math::RealFunction& potential, const Math::RealFunction& potentialDeriv, double phi0, double phiDot0, double phiEnd, double a0, int curvature = 0);

    const Math::TableFunction<double, double>& phi() const { return phi_; }
    const Math::TableFunction<double, double>& phiDot() const { return phiDot_; }
    const Math::TableFunction<double, double>& H() const { return H_; }
    const Math::TableFunction<double, double>& a() const { return a_; }
    const Math::TableFunction<double, double>& logA() const { return logA_; }
    const Math::TableFunction<double, double>& z() const { return z_; }
    const Math::TableFunction<double, double>& zPrime() const { return zPrime_; }
    const Math::TableFunction<double, double>& zDoublePrime() const { return zDoublePrime_; }
    const Math::TableFunction<double, double>& tau() const { return tau_; }

    double scalarP(double k) const;

private:
    Math::TableFunction<double, double> phi_, phiDot_, H_, a_, logA_, z_, zPrime_, zDoublePrime_, tau_, zRatioTau_, aHInverse_;
    int K_;
    double tauStart_, tStart_;
    double aHStart_, aHEnd_;
};

#endif

