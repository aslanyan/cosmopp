#ifndef COSMO_PP_INFLATION_HPP
#define COSMO_PP_INFLATION_HPP

class SingleFieldInflation
{
public:
    SingleFieldInflation(const Math::RealFunction& potential, const Math::RealFunction& potentialDeriv, double phi0, double phiDot0, double a0, double aDot0, double curvature = 0);
};

#endif

