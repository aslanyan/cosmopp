#ifndef COSMO_CPP_PARAMETRIC_FUNCTION_HPP
#define COSMO_CPP_PARAMETRIC_FUNCTION_HPP

#include <vector>

#include <function.hpp>
#include <macros.hpp>

namespace Math
{

/// A parametric function class.

/// This is an abstract class for parametric functions. Can be used for curve fitting, for example. The templete argument is the number of the parameters.
template<int N>
class ParametricFunction : public Function<double, double>
{
public:
    /// Constructor.

    /// Constructor, all the parameters are given value 0 by default.
    ParametricFunction() : params_(N, 0) {}

    /// Destructor.
    virtual ~ParametricFunction() {}
    
    /// Parameter access function.
    /// \param i The index of the parameter (from 0 to N - 1).
    /// \return The value of that parameter.
    virtual double parameter(int i) const { check(i >= 0 && i < N, "invalid i = " << i); return params_[i]; }

    /// Parameter access function.
    /// \param i The index of the parameter (from 0 to N - 1).
    /// \return A reference to the parameter (can be used to change it).
    virtual double& parameter(int i) { check(i >= 0 && i < N, "invalid i = " << i); return params_[i]; }
    
    /// A purely virtual function to evaluate the function.
    /// \param x The argument of the function.
    /// \return The value of the function.
    virtual double evaluate(double x) const = 0;
protected:
    std::vector<double> params_;
};
    
} //namespace Math

#endif
