#ifndef COSMO_PP_CHI_SQUARED_HPP
#define COSMO_PP_CHI_SQUARED_HPP

#include <cmath>

#include <macros.hpp>
#include <function.hpp>

namespace Math
{

class ChiSquared : public RealFunction
{
public:
    ChiSquared(int dof)
    {
        check(dof >= 1, "invalid number of degrees of freedom " << dof << ", needs to be positive");
        factor_ = double(dof) / 2 - 1;
        add_ = -std::lgamma(double(dof) / 2) - double(dof) / 2 * std::log(2.0);
    }

    double evaluate(double x) const
    {
        check(x >= 0, "invalid value of x = " << x << ", needs to be non-negative");

        if(x == 0)
            return 0;

        double res = add_ + factor_ * std::log(x) - x / 2.0;
        return std::exp(res);
    }

private:
    double factor_;
    double add_;
};

} // namespace Math

#endif

