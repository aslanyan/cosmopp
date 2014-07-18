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

} // namespace Math

#endif

