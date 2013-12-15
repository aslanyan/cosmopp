#ifndef COSMO_PP_POLYNOMIAL_HPP
#define COSMO_PP_POLYNOMIAL_HPP

#include <parametric_function.hpp>

namespace Math
{

/// A polynomial of degree N - 1 (it has N parameters as coefficients).
template<int N>
class Polynomial : public ParametricFunction<N>
{
public:
    /// Constructor.

    /// Constructor, all of the coefficients are set to 0 by default.
    Polynomial() : ParametricFunction<N>() {}

    /// Destructor.
    virtual ~Polynomial() {}
    
    /// Evaluate the polynomial.
    /// \param x The argument.
    /// \return The value of the polynomial.
    virtual double evaluate(double x) const
    {
        double res = 0;
        double current = 1;
        for(int i = 0; i < N; ++i)
        {
            res += ParametricFunction<N>::parameter(i) * current;
            current *= x;
        }
        
        return res;
    }
};
    
} //namespace Math

#endif
