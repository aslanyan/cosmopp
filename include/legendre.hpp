#ifndef COSMO_PP_LEGENDRE_HPP
#define COSMO_PP_LEGENDRE_HPP

#include <vector>
#include <string>
#include <sstream>
#include <cmath>

namespace Math
{

class Legendre
{
public:
    Legendre() : vals_(10000, 1) {}
    ~Legendre() {}

    double calculate(unsigned int l, double x)
    {
        if(vals_.size() < l + 1)
            vals_.resize(l + 1);

        //vals_[0] = 1;
        vals_[1] = x;
        for(int l1 = 2; l1 <= l; ++l1)
            vals_[l1] = (2 - 1.0 / l1) * x * vals_[l1 - 1] - (1 - 1.0 / l1) * vals_[l1 - 2];

        return vals_[l];
    }

private:
    std::vector<double> vals_;
};

class AssociatedLegendre
{
public:
    AssociatedLegendre() : vals_(10000, 1) {}
    ~AssociatedLegendre() {}

    double calculate(unsigned int l, int m, double x)
    {
        if(int(l) < m || int(l) < -m)
            return 0;

        if(l == 0)
            return 1;

        if(m < 0)
        {
            m = -m;
            double factor = (m % 2 ? -1.0 : 1.0);
            for(int i = l - m + 1; i <= l + m; ++i)
                factor /= i;

            return factor * calculate(l, m, x);
        }

        if(vals_.size() < l + 1)
            vals_.resize(l + 1);

        if(int(l) == m)
        {

            const double f = std::sqrt(1.0 - x * x);

            vals_[0] = 1.0;

            for(int l1 = 1; l1 <= l; ++l1)
                vals_[l1] = -(2 * l1 - 1) * f * vals_[l1 - 1];

            return vals_[l];
        }

        vals_[m] = calculate(m, m, x);
        vals_[m + 1] = (2 * m + 1) * x * vals_[m];
        
        for(int l1 = m + 2; l1 <= l; ++l1)
            vals_[l1] = ((2 * l1 - 1) * x * vals_[l1 - 1] - (l1 + m - 1) * vals_[l1 - 2]) / (l1 - m);
        
        return vals_[l];
    }

private:
    std::vector<double> vals_;
};

} // namespace Math

#endif

