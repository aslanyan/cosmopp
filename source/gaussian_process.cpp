#include <macros.hpp>
#include <gaussian_process.hpp>

#include <gmd.h>
#include <lavd.h>
#include <laslv.h>
#include <lavli.h>
#include <blas2pp.h>
#include <blas3pp.h>

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
GaussianProcess::set(const std::vector<std::vector<double> >& x, const std::vector<double>& y)
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

    LUFactorizeIP(k_, pivotK_);
    LaLUInverseIP(k_, pivotK_);
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
    Blas_Mat_Mat_Mult(prod1, prod1, prod2, false, true);

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
