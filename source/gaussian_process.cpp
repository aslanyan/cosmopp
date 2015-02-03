#include <fstream>
#include <string>
#include <cmath>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <numerics.hpp>
#include <timer.hpp>
#include <gaussian_process.hpp>

#include <gmd.h>
#include <lavd.h>
#include <laslv.h>
#include <lavli.h>
#include <blas2pp.h>
#include <blas3pp.h>

#include "Minuit2/FCNGradientBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"

namespace
{

class GaussianProcessFCN : public ROOT::Minuit2::FCNBase
{
public:
    GaussianProcessFCN(const std::vector<std::vector<double> >& x, const std::vector<double>& y) : x_(x), y_(y)
    {
        check(x_.size() == y_.size(), "");
        check(!x_.empty(), "");
    }

    virtual double operator()(const std::vector<double>& par) const
    {
        LaGenMatDouble k(x_.size(), x_.size());
        LaGenMatDouble kInv(x_.size(), x_.size());
        double logDet;

        calculateCovarianceInverse(par, k, kInv, logDet);

        double res = logDet;

        for(int i = 0; i < x_.size(); ++i)
        {
            for(int j = 0; j < x_.size(); ++j)
                res += kInv(i, j) * y_[i] * y_[j];
        }

        return res;
    }

    virtual std::vector<double> Gradient(const std::vector<double>& par) const
    {
        LaGenMatDouble k(x_.size(), x_.size());
        LaGenMatDouble kInv(x_.size(), x_.size());
        double logDet;

        calculateCovarianceInverse(par, k, kInv, logDet);

        check(!par.empty(), "");
        std::vector<double> grad(par.size());

        std::vector<double> alpha(y_.size(), 0);
        for(int i = 0; i < y_.size(); ++i)
        {
            for(int j = 0; j < y_.size(); ++j)
                alpha[i] += kInv(i, j) * y_[j];
        }

        for(int i = 0; i < grad.size(); ++i)
        {
            LaGenMatDouble kDeriv = k;
            for(int j = 0; j < x_.size(); ++j)
            {
                for(int l = 0; l < x_.size(); ++l)
                {
                    double factor;
                    if(i == grad.size() - 1)
                        factor = 2.0 / par[i];
                    else
                    {
                        const double delta = x_[j][i] - x_[l][i];
                        factor = delta * delta / (par[i] * par[i] * par[i]);
                    }

                    kDeriv(i, j) *= factor;
                }
            }

            double trace = 0;
            for(int j = 0; j < x_.size(); ++j)
            {
                for(int l = 0; l < x_.size(); ++l)
                {
                    trace += alpha[j] * alpha[l] * kDeriv(l, j);
                    trace -= kInv(j, l) * kDeriv(l, j);
                }
            }
            grad[i] = -trace;

#ifdef CHECKS_ON
            const double epsilon = 1e-5;
            std::vector<double> par1 = par;
            par1[i] += epsilon;
            const double gradCheck = ((*this)(par1) - (*this)(par)) / epsilon;
            check(Math::areEqual(gradCheck, grad[i], 1e-1), "gradient discrepancy for parameter " << i << ", expected " << gradCheck << " obtained " << grad[i]);
#endif
        }

        return grad;
    }
    
    double Up() const {return 1.;}
    
private:
    double cov(const std::vector<double>& x1, const std::vector<double>& x2, double sigma2, const std::vector<double>& l2) const
    {
        check(!l2.empty(), "");
        check(x1.size() == l2.size(), "");
        check(x2.size() == l2.size(), "");

        check(sigma2 > 0, "");

        double d = 0;

        for(int i = 0; i < l2.size(); ++i)
        {
            check(l2[i] > 0, "");
            const double delta = x1[i] - x2[i];

            d += delta * delta / l2[i];
        }

        return sigma2 * std::exp(- d / 2);
    }

    void calculateCovarianceInverse(const std::vector<double>& par, LaGenMatDouble& k, LaGenMatDouble& kInv, double& logDet) const
    {
        check(x_.size() == y_.size(), "");
        check(!x_.empty(), "");

        const int n = x_[0].size();

        check(par.size() == n + 1, "");

        double sigma2 = par[n];

        check(sigma2 > 0, "");
        sigma2 = sigma2 * sigma2;

        std::vector<double> l2 = par;
        l2.pop_back();

        for(int i = 0; i < n; ++i)
        {
            check(l2[i] > 0, "");
            l2[i] = l2[i] * l2[i];
        }


        for(int i = 0; i < x_.size(); ++i)
        {
            for(int j = i; j < x_.size(); ++j)
            {
                const double c = cov(x_[i], x_[j], sigma2, l2);
                k(i, j) = c;
                k(j, i) = c;
            }
        }

        //hack
        for(int i = 0; i < x_.size(); ++i)
            k(i, i) *= 1.00001;

        /*
        std::ofstream out("cov_matrix.txt");
        for(int i = 0; i < x_.size(); ++i)
        {
            for(int j = 0; j < x_.size(); ++j)
                out << k(i, j) << '\t';
            out << std::endl;
        }
        out.close();
        */

        kInv = k;

        LaVectorLongInt pivotK(x_.size());

        LUFactorizeIP(kInv, pivotK);

        logDet = 0;
        int signC = 1;

        StandardException exc;
        for(int i = 0; i < x_.size(); ++i)
        {
            if(Math::areEqual(kInv(i, i), 0.0, 1e-15))
            {
                std::string exceptionStr = "The determinant of the covariance matrix is 0. The covariance matrix must be positive definite.";
                exc.set(exceptionStr);
                throw exc;
            }

            signC *= (kInv(i, i) < 0 ? -1 : 1);
            logDet += std::log(std::abs(kInv(i, i)));
        }

        //check(signC == 1, "the covariance matrix determinant is negative!");

        LaLUInverseIP(kInv, pivotK);

#ifdef CHECKS_ON
        const double epsilon = 1e-5;
        LaGenMatDouble prod(x_.size(), x_.size());
        Blas_Mat_Mat_Mult(k, kInv, prod);

        for(int i = 0; i < x_.size(); ++i)
        {
            for(int j = 0; j < x_.size(); ++j)
            {
                if(i == j)
                {
                    check(Math::areEqual(prod(i, i), 1.0, epsilon), "diagonal element (" << i << ", " << i << ") should be 1 but it is " << prod(i, i));
                }
                else
                    check(Math::areEqual(prod(i, j), 0.0, epsilon), "non-diagonal element (" << i << ", " << j << ") should be 0 but it is " << prod(i, j));
            }
        }
#endif
    }

private:
    const std::vector<std::vector<double> >& x_;
    const std::vector<double>& y_;
};

} // namespace

namespace Math
{

GaussianProcess::GaussianProcess(int dim, CovarianceType covType) : dim_(dim), covType_(covType)
{
    check(dim_ > 0, "invalid dimension " << dim_ << ", should be positive");
    check(covType_ >= SQUARED_EXP && covType_ < COVARIANCE_TYPE_MAX, "invalid covariance type " << covType_);

    switch(covType)
    {
    case SQUARED_EXP:
        sigma2_ = 1;
        l2_.resize(dim_, 1);
        break;

    default:
        check(false, "");
    }
}

void
GaussianProcess::setParams(const std::vector<double>& params)
{
    switch(covType_)
    {
    case SQUARED_EXP:
        check(params.size() >= dim_ + 1, "need at least" << dim_ + 1 << " params");
        check(l2_.size() == dim_, "");

        for(int i = 0; i < dim_; ++i)
        {
            check(params[i] > 0, "");
            l2_[i] = params[i] * params[i];
        }

        check(params[dim_] > 0, "");
        sigma2_ = params[dim_] * params[dim_];

        break;
    default:
        check(false, "");
    }
}

void
GaussianProcess::set(const std::vector<std::vector<double> >& x, const std::vector<double>& y, bool calculateParamsByMaxLike)
{
    check(!x.empty(), "x is empty");
    check(x.size() == y.size(), "the sizes of x and y don't match");

    y_ = y;

    x_.resize(x.size());
    for(int i = 0; i < x_.size(); ++i)
    {
        check(x[i].size() == dim_, "the size of input " << i << " is " << x[i].size() << " so it doesn't match the dimension " << dim_);
        x_[i] = x[i];
    }

    if(calculateParamsByMaxLike)
    {
        std::vector<std::vector<double> > x1 = x_;
        std::vector<double> y1 = y_;

        // setting max size
        const int maxSize = 1000;
        if(x1.size() > maxSize)
        {
            x1.resize(maxSize);
            y1.resize(maxSize);
        }

        GaussianProcessFCN fcn(x1, y1);
        std::vector<double> min(dim_ + 1, 1e-5);
        std::vector<double> max(dim_ + 1, 10);
        max[dim_] = 1e3;
        std::vector<double> error(dim_ + 1, 1e-3);

        std::vector<double> starting(dim_ + 1);
        for(int i = 0; i < dim_; ++i)
        {
            starting[i] = std::sqrt(l2_[i]);
            if(starting[i] < min[i])
                starting[i] = min[i];
            if(starting[i] > max[i])
                starting[i] = max[i];
        }

        starting[dim_] = std::sqrt(sigma2_);
        if(starting[dim_] < min[dim_])
            starting[dim_] = min[dim_];
        if(starting[dim_] > max[dim_])
            starting[dim_] = max[dim_];

        ROOT::Minuit2::MnUserParameters upar;
        for(int i = 0; i <= dim_; ++i)
        {
            std::stringstream paramName;
            if(i < dim_)
                paramName << "l_" << i;
            else
                paramName << "sigma";

            upar.Add(paramName.str().c_str(), starting[i], error[i], min[i], max[i]);
        }

        Timer t1("GAUSSIAN PROCESS MINIMIZATION");
        t1.start();
        ROOT::Minuit2::MnMigrad migrad(fcn, upar);
        ROOT::Minuit2::FunctionMinimum minRes = migrad();
        ROOT::Minuit2::MnUserParameters result = minRes.UserParameters();
        t1.end();
        
        for(int i = 0; i <= dim_; ++i)
        {
            std::stringstream paramName;
            if(i < dim_)
                paramName << "l_" << i;
            else
                paramName << "sigma";

            const double param = result.Value(paramName.str().c_str());
            const double e = result.Error(paramName.str().c_str());

            output_screen(paramName.str() << " = " << param << " +/- " << e << std::endl);

            if(i < dim_)
                l2_[i] = param * param;
            else
                sigma2_ = param * param;
        }
    }

    k_.resize(x_.size(), x_.size());

    for(int i = 0; i < x_.size(); ++i)
    {
        for(int j = i; j < x_.size(); ++j)
        {
            const double c = cov(x_[i], x_[j]);
            k_(i, j) = c;
            k_(j, i) = c;
        }
    }

    //hack
    for(int i = 0; i < x_.size(); ++i)
        k_(i, i) *= 1.00001;

    LaGenMatDouble kOrig = k_;

    LaVectorLongInt pivotK(x_.size());

    LUFactorizeIP(k_, pivotK);
    LaLUInverseIP(k_, pivotK);

#ifdef CHECKS_ON
    const double epsilon = 1e-5;
    LaGenMatDouble prod(x_.size(), x_.size());
    Blas_Mat_Mat_Mult(kOrig, k_, prod);

    for(int i = 0; i < x_.size(); ++i)
    {
        for(int j = 0; j < x_.size(); ++j)
        {
            if(i == j)
            {
                check(Math::areEqual(prod(i, i), 1.0, epsilon), "diagonal element (" << i << ", " << i << ") should be 1 but it is " << prod(i, i));
            }
            else
                check(Math::areEqual(prod(i, j), 0.0, epsilon), "non-diagonal element (" << i << ", " << j << ") should be 0 but it is " << prod(i, j));
        }
    }
#endif
}

void
GaussianProcess::calculate(const std::vector<std::vector<double> >& input, std::vector<double>& mean, LaGenMatDouble& covariance) const
{
    check(!input.empty(), "empty input given");

    const int n = input.size();
    const int m = x_.size();

    check(m > 0, "");
    check(y_.size() == m, "");

    LaGenMatDouble c(n, m);

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < m; ++j)
            c(i, j) = cov(input[i], x_[j]);
    }
    
    LaGenMatDouble prod1(n, m);
    Blas_Mat_Mat_Mult(c, k_, prod1, false);

    mean.resize(n);
    for(int i = 0; i < n; ++i)
    {
        mean[i] = 0;
        for(int j = 0; j < m; ++j)
            mean[i] += prod1(i, j) * y_[j];
    }

    LaGenMatDouble prod2(n, n);
    Blas_Mat_Mat_Mult(prod1, c, prod2, false, true);

    covariance.resize(n, n);
    for(int i = 0; i < n; ++i)
    {
        for(int j = i; j < n; ++j)
        {
            const double c = cov(input[i], input[j]) - prod2(i, j);
            covariance(i, j) = c;
            covariance(j, i) = c;
        }
    }
}

double
GaussianProcess::cov(const std::vector<double>& x1, const std::vector<double>& x2) const
{
    check(x1.size() == dim_, "");
    check(x2.size() == dim_, "");

    check(l2_.size() == dim_, "");
    check(sigma2_ > 0, "");

    double d = 0;

    for(int i = 0; i < dim_; ++i)
    {
        check(l2_[i] > 0, "");
        const double delta = x1[i] - x2[i];

        d += delta * delta / l2_[i];
    }

    return sigma2_ * std::exp(- d / 2);
}

} // namespace Math

