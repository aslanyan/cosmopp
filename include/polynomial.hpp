#ifndef COSMO_PP_POLYNOMIAL_HPP
#define COSMO_PP_POLYNOMIAL_HPP

#include <vector>

#include <parametric_function.hpp>

namespace Math
{

/// A polynomial of degree n - 1 (it has n parameters as coefficients).
class Polynomial : public ParametricFunction
{
public:
    /// Constructor.

    /// Constructor, all of the coefficients are set to 0 by default.
    Polynomial(unsigned int n) : ParametricFunction(n) {}

    /// Constructor.
    
    /// Copy constructor.
    /// \param other The polynomial to copy from.
    Polynomial(const Polynomial& other) : ParametricFunction(other) {}

    /// Destructor.
    virtual ~Polynomial() {}

    Polynomial& operator = (const Polynomial& other)
    {
        ParametricFunction::operator = (other);
        return *this;
    }
    
    /// Evaluate the polynomial.
    /// \param x The argument.
    /// \return The value of the polynomial.
    virtual double evaluate(double x) const
    {
        double res = 0;
        double current = 1;
        for(int i = 0; i < numberOfParams(); ++i)
        {
            res += ParametricFunction::parameter(i) * current;
            current *= x;
        }
        
        return res;
    }

    Polynomial operator + (const Polynomial& other) const
    {
        unsigned int n = numberOfParams();
        if(other.numberOfParams() > n)
            n = other.numberOfParams();

        std::vector<double> coefficients(n, 0.0);

        for(int i = 0; i < numberOfParams(); ++i)
            coefficients[i] += parameter(i);
        for(int i = 0; i < other.numberOfParams(); ++i)
            coefficients[i] += other.parameter(i);
        
        clearZeros(coefficients);

        Polynomial result(coefficients.size());
        for(int i = 0; i < coefficients.size(); ++i)
            result.parameter(i) = coefficients[i];

        return result;
    }

    Polynomial operator - (const Polynomial& other) const
    {
        unsigned int n = numberOfParams();
        if(other.numberOfParams() > n)
            n = other.numberOfParams();

        std::vector<double> coefficients(n, 0.0);

        for(int i = 0; i < numberOfParams(); ++i)
            coefficients[i] += parameter(i);
        for(int i = 0; i < other.numberOfParams(); ++i)
            coefficients[i] -= other.parameter(i);
        
        clearZeros(coefficients);

        Polynomial result(coefficients.size());
        for(int i = 0; i < coefficients.size(); ++i)
            result.parameter(i) = coefficients[i];

        return result;
    }

    Polynomial operator * (const Polynomial& other) const
    {
        unsigned int n = numberOfParams() + other.numberOfParams() - 1;
        std::vector<double> coefficients(n, 0.0);

        for(int i = 0; i < numberOfParams(); ++i)
        {
            for(int j = 0; j < other.numberOfParams(); ++j)
            {
                int k = i + j;
                check(k < n, "");
                coefficients[k] += parameter(i) * other.parameter(j);
            }
        }

        clearZeros(coefficients);

        Polynomial result(coefficients.size());
        for(int i = 0; i < coefficients.size(); ++i)
            result.parameter(i) = coefficients[i];

        return result;
    }
    
    Polynomial operator * (double c) const
    {
        Polynomial res(numberOfParams());
        for(int i = 0; i < numberOfParams(); ++i)
            res.parameter(i) = c * parameter(i);

        return res;
    }

    Polynomial operator / (double c) const
    {
        check(c != 0, "cannot divide by 0");

        Polynomial res(numberOfParams());
        for(int i = 0; i < numberOfParams(); ++i)
            res.parameter(i) = parameter(i) / c;

        return res;
    }

private:
    void clearZeros(std::vector<double>& v) const
    {
        while(!v.empty())
        {
            int i = v.size() - 1;
            if(v[i] != 0.0)
                break;
            v.pop_back();
        }
    }
};

inline
Polynomial operator *(double c, const Polynomial& p)
{
    Polynomial res = p;
    for(int i = 0; i < res.numberOfParams(); ++i)
        res.parameter(i) *= c;

    return res;
}
    
} //namespace Math

#endif
