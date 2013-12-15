#ifndef COSMO_PP_SPIN_SPHERICAL_HARMONICS_HPP
#define COSMO_PP_SPIN_SPHERICAL_HARMONICS_HPP

#include <cmath>

#include <macros.hpp>
#include <angular_coordinates.hpp>
#include <complex_types.hpp>
#include <math_constants.hpp>

namespace Math
{
    
/// A function to calculate the spin spherical harmonic.
/// \param s The spin.
/// \param l l-value.
/// \param m m-value.
/// \param angularCoords The angular coordinates.
/// \return The spin spherical harmonic value.
inline ComplexDouble spinSphericalHarmonic(int s, int l, int m, AngularCoordinates angularCoords)
{
    check(l >= 0, "invalid l");
    check(l >= m && l >= -m, "invalid m");
    check(l >= s && l >= -s, "invalid s");
    
    //calculating (l + m)! / (l + s)!
    double factorialRatio1 = 1;
    if(l + m > l + s)
    {
        for(int i = l + s + 1; i <= l + m; ++i)
            factorialRatio1 *= i;
    }
    else
    {
        for(int i = l + m + 1; i <= l + s; ++i)
            factorialRatio1 /= i;
    }
    
    //calculating (l - m)! / (l - s)!
    double factorialRatio2 = 1;
    if(l - m > l - s)
    {
        for(int i = l - s + 1; i <= l - m; ++i)
            factorialRatio2 *= i;
    }
    else
    {
        for(int i = l - m + 1; i <= l - s; ++i)
            factorialRatio2 /= i;
    }
    
    double factor = factorialRatio1 * factorialRatio2 * (2 * l + 1) / (4 * pi);
    factor = std::sqrt(factor);
    
    const double sinThetaHalf = std::sin(angularCoords.theta / 2);
    const double cosThetaHalf = std::cos(angularCoords.theta / 2);
    
    const double sinMPhi = std::sin(m * angularCoords.phi);
    const double cosMPhi = std::cos(m * angularCoords.phi);
    
    int rMin = 0;
    if(m - s > 0)
        rMin = m - s;
    
    int rMax = l - s;
    if(l + m < rMax)
        rMax = l + m;
    check(rMax >= rMin, "");
    
    //f0 = (-1)^(l + m - r - s)
    int f0 = 1;
    if((l + m - rMin - s) % 2)
        f0 = -1;
    
    //f1 = binomial(l - s, r)
    double f1 = 1;
    for(int i = 1; i <= l - s - rMin; ++i)
    {
        f1 *= (rMin + i);
        f1 /= i;
    }
    
    //f2 = binomial(l + s, r + s - m)
    double f2 = 1;
    const int a = l + s, b = rMin + s - m;
    check(a >= 0 && b >= 0 , "");
    check(a >= b, "");
    for(int i = 1; i <= a - b; ++i)
    {
        f2 *= (b + i);
        f2 /= i;
    }
    
    //f3 = (cos(theta / 2))^(2r + s - m)
    check(2 * rMin + s - m >= 0, "");
    double f3 = 1;
    for(int i = 1; i <= 2 * rMin + s - m; ++i)
        f3 *= cosThetaHalf;
    
    double sum = 0;
    for(int r = rMin; r <= rMax; ++r)
    {
        //f4 = ((sin(theta / s))^(2l - 2r - s + m)
        double f4 = 1;
        check(2 * l - 2 * r - s + m >= 0, "");
        for(int i = 1; i <= 2 * l - 2 * r - s + m; ++i)
            f4 *= sinThetaHalf;
        
        sum += f0 * f1 * f2 * f3 * f4;
        
        if(r == rMax)
            break;
        
        //prepare f0 -> f3 for next step
        f0 *= -1;
        
        f1 *= (l - s - r);
        f1 /= (r + 1);
        
        f2 *= (l + m - r);
        f2 /= (r + 1 + s - m);
        
        f3 *= (cosThetaHalf * cosThetaHalf);
    }
    
    ComplexDouble result(factor * sum * cosMPhi, factor * sum * sinMPhi);
    return result;
}
	
} //namespace Math

#endif
