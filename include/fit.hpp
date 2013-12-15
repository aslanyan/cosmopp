#ifndef COSMO_PP_FIT_HPP
#define COSMO_PP_FIT_HPP

#include <vector>
#include <sstream>
#include <string>

#include <parametric_function.hpp>
#include <macros.hpp>

#include "Minuit2/FCNGradientBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"

namespace Math
{

template<int N>
class Chi2Calculator : public ROOT::Minuit2::FCNBase
{
public:
    Chi2Calculator(ParametricFunction<N>& f, const std::vector<double>& x, const::std::vector<double>& y) : f_(f), x_(x), y_(y)
    {
        check(x_.size() == y_.size(), "");
        check(!x_.empty(), "");
    }
    
    virtual double operator()(const std::vector<double>& par) const
    {
        check(par.size() == N, "");
        for(int i = 0; i < N; ++i)
            f_.parameter(i) = par[i];
        
        double res = 0;
        for(int i = 0; i < x_.size(); ++i)
        {
            const double diff = y_[i] - f_.evaluate(x_[i]);
            res += diff * diff;
        }
        return res;
    }
    
    double Up() const {return 1.;}
    
private:
    ParametricFunction<N>& f_;
    std::vector<double> x_;
    std::vector<double> y_;
};

/// Curve fitting class.

/// This class provides the functionality for curve fitting. Minimization is performed using the Minuit C++ library, this needs to be installed.
template<int N>
class Fit
{
public:
    /// Constructor.

    /// Constructs the curve fitter.
    /// \param f The parametric function used for fitting the data.
    /// \param x A vector containing the x coordinates of the data points.
    /// \param y A vector containing the y coordinates of the data points. Must have the same size as x.
    /// \param starting A vector containing the starting values of the parameters.
    /// \param error A vector containing the errors of the parameters.
    /// \param min A vector containing the lower limits for the ranges of the parameters. Provide an empty vector to allow unlimited ranges.
    /// \param max A vector containing the upper limits for the ranges of the parameters.
    Fit(ParametricFunction<N>& f, const std::vector<double>& x, const::std::vector<double>& y, const std::vector<double>& starting, const std::vector<double>& error, const std::vector<double>& min, const std::vector<double>& max) : f_(f)
    {
        calc_ = new Chi2Calculator<N>(f, x, y);
        check(starting.size() == N, "");
        check(error.size() == N, "");
        check(min.size() == N || min.empty(), "");
        check(max.size() == N || max.empty(), "");
        check(min.size() == max.size(), "");
        
        for(int i = 0; i < N; ++i)
        {
            std::stringstream paramName;
            paramName << "param_" << i;
            if(min.empty())
                upar_.Add(paramName.str().c_str(), starting[i], error[i]);
            else
                upar_.Add(paramName.str().c_str(), starting[i], error[i], min[i], max[i]);
        }
    }
    
    /// Destructor.
    ~Fit() { delete calc_; }
    
    /// Does the actual fitting.
    /// \param params A vector that will contain the best fit values of the parameters upon return.
    /// \param error A vector that will contain the uncertainties of the values in params upon return.
    void fit(std::vector<double>& params, std::vector<double>& error)
    {
        ROOT::Minuit2::MnMigrad migrad(*calc_, upar_);
        ROOT::Minuit2::FunctionMinimum minRes = migrad();
        ROOT::Minuit2::MnUserParameters result = minRes.UserParameters();
        
        params.clear();
        error.clear();
        
        for(int i = 0; i < N; ++i)
        {
            std::stringstream paramName;
            paramName << "param_" << i;
            const double param = result.Value(paramName.str().c_str());
            f_.parameter(i) = param;
            params.push_back(param);
            error.push_back(result.Error(paramName.str().c_str()));
        }
    }
    
private:
    Chi2Calculator<N>* calc_;
    ROOT::Minuit2::MnUserParameters upar_;
    ParametricFunction<N>& f_;
};
    
} //namespace Math

#endif
