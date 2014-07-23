#ifndef COSMO_PP_SPHERICAL_HARMONICS_HPP
#define COSMO_PP_SPHERICAL_HARMONICS_HPP

#include <cmath>

#include <macros.hpp>
#include <complex_types.hpp>
#include <legendre.hpp>
#include <math_constants.hpp>

namespace Math
{

class SphericalHarmonics
{
public:
    SphericalHarmonics() {}
    ~SphericalHarmonics() {}

    ComplexDouble calculate(unsigned int l, int m, double theta, double phi)
    {
        check(int(l) >= m && int(l) >= -m, "");

        const double ctheta = std::cos(theta);
        const double associatedLeg = leg_.calculate(l, m, ctheta);

        double f = double(2 * l + 1) / (4 * pi);

        if(m >= 0)
        {
            for(int i = int(l) - m + 1; i  <= int(l) + m; ++i)
                f /= double(i);
        }

        else
        {
            for(int i = int(l) + m + 1; i <= int(l) - m; ++i)
                f *= double(i);
        }

        f = std::sqrt(f) * associatedLeg;
        ComplexDouble res(f * std::cos(m * phi), f * std::sin(m * phi));
        return res;
    }

private:
    AssociatedLegendre leg_;
};

} // namespace Math
#endif

