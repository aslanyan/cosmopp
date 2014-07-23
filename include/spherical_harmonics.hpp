#ifndef COSMO_PP_SPHERICAL_HARMONICS_HPP
#define COSMO_PP_SPHERICAL_HARMONICS_HPP

#include <vector>
#include <cmath>

#include <complex_types.hpp>
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
        if(int(l) < m || int(l) < -m)
            return 0;

        const double f = std::sqrt(double(2 * l + 1) / (4 * pi)) * modifiedAssociatedLeg(l, m, std::cos(theta));
        ComplexDouble res(f * std::cos(m * phi), f * std::sin(m * phi));
        return res;
    }

private:
    double modifiedAssociatedLeg(unsigned int l, int m, double x)
    {
        if(int(l) < m || int(l) < -m)
            return 0;

        if(l == 0)
            return 1;

        if(m < 0)
        {
            m = -m;
            double factor = (m % 2 ? -1.0 : 1.0);
            return factor * modifiedAssociatedLeg(l, m, x);
        }

        if(vals_.size() < l + 1)
            vals_.resize(l + 1);

        if(int(l) == m)
        {

            const double f = 1.0 - x * x;

            vals_[0] = 1.0;

            for(int l1 = 1; l1 <= l; ++l1)
                vals_[l1] = -std::sqrt(f * (1.0 - 1.0 / double(2 * l1))) * vals_[l1 - 1];

            return vals_[l];
        }

        vals_[m] = modifiedAssociatedLeg(m, m, x);
        vals_[m + 1] = std::sqrt(2 * m + 1) * x * vals_[m];
        
        for(int l1 = m + 2; l1 <= l; ++l1)
            vals_[l1] = std::sqrt(double(4 * l1 * l1 - 4 * l1 + 1) / double(l1 * l1 - m * m)) * x * vals_[l1 - 1] - std::sqrt((1.0 - 1.0 / double(l1 + m)) * (1.0 - 1.0 / double(l1 - m))) * vals_[l1 - 2];
        
        return vals_[l];
    }
private:
    std::vector<double> vals_;
};

} // namespace Math
#endif

