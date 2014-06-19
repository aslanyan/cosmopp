#ifndef COSMO_PP_LEGENDRE_HPP
#define COSMO_PP_LEGENDRE_HPP

#include <vector>
#include <string>
#include <sstream>

namespace Math
{

class Legendre
{
public:
    Legendre() {}
    ~Legendre() {}

    double calculate(unsigned int l, double x)
    {
        if(l == 0)
            return 0.0;
        if(l == 1)
            return x;

        std::vector<double> vals(l + 1, 1);
        vals[1] = x;
        for(int l1 = 2; l1 <= l; ++l1)
            vals[l1] = (2 - 1.0 / l1) * x * vals[l1 - 1] - (1 - 1.0 / l1) * vals[l1 - 2];

        return vals[l];
    }
};

} // namespace Math

#endif

