#ifndef COSMO_PP_PARAMETRIC_FUNCTION_HPP
#define COSMO_PP_PARAMETRIC_FUNCTION_HPP

#include <vector>

#include <function.hpp>
#include <macros.hpp>

namespace Math
{

/// A parametric function class.

/// This is an abstract class for parametric functions. Can be used for curve fitting, for example.
class ParametricFunction : public Function<double, double>
{
public:
    /// Constructor.

    /// Constructor, all the parameters are given value 0 by default.
    /// \param n The number of parameters
    ParametricFunction(unsigned int n) : params_(n, 0) {}

    /// Constructor.

    /// Copy constructor.
    /// \param other The object to copy from.
    ParametricFunction(const ParametricFunction& other) : params_(other.params_) {}

    /// Destructor.
    virtual ~ParametricFunction() {}

    /// Copy operator from another object.
    /// \param other Object to be copied from.
    /// \return Reference to self after the copy.
    ParametricFunction& operator = (const ParametricFunction& other)
    {
        params_ = other.params_;
        return *this;
    }

    /// Number of parameters.
    /// \return The number of parameters.
    virtual double numberOfParams() const { return params_.size(); }
    
    /// Parameter access function.
    /// \param i The index of the parameter (from 0 to n - 1).
    /// \return The value of that parameter.
    virtual double parameter(int i) const { check(i >= 0 && i < params_.size(), "invalid i = " << i); return params_[i]; }

    /// Parameter access function.
    /// \param i The index of the parameter (from 0 to n - 1).
    /// \return A reference to the parameter (can be used to change it).
    virtual double& parameter(int i) { check(i >= 0 && i < params_.size(), "invalid i = " << i); return params_[i]; }
    
    /// A purely virtual function to evaluate the function.
    /// \param x The argument of the function.
    /// \return The value of the function.
    virtual double evaluate(double x) const = 0;

protected:
    std::vector<double> params_;
};
    
} //namespace Math

#endif

