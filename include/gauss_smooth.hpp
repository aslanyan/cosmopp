#ifndef COSMO_PP_GAUSS_SMOOTH
#define COSMO_PP_GAUSS_SMOOTH

#include <vector>
#include <algorithm>
#include <cmath>

#include <macros.hpp>
#include <function.hpp>

namespace Math
{

class GaussSmooth : public RealFunction
{
public:
    inline GaussSmooth(const std::vector<double>& x, const std::vector<double>& y, double sigma);

    ~GaussSmooth() {}

    inline virtual double evaluate(double x) const;

private:
    inline double kernel(double x, double x1) const;

private:
    const std::vector<double> x_, y_;
    const double sigmaSq_;
    const double sigma_;
};

GaussSmooth::GaussSmooth(const std::vector<double>& x, const std::vector<double>& y, double sigma) : x_(x), y_(y), sigmaSq_(sigma * sigma), sigma_(sigma)
{
    check(!x_.empty(), "need to have at least 1 point");
    check(x_.size() == y_.size(), "the sizes of vectors x and y must match");
    check(sigma > 0, "invalid sigma = " << sigma << ", must be positive");

#ifdef CHECKS_ON
    //checking that the x vector is sorted
    for(unsigned long i = 1; i < x_.size(); ++i)
    {
        check(x_[i] >= x_[i - 1], "the x vector is not sorted");
    }
#endif
}

double
GaussSmooth::evaluate(double x) const
{
    check(x_.size() == y_.size(), "");
    double res = 0, norm = 0;
    const std::vector<double>::const_iterator lower = std::lower_bound(x_.begin(), x_.end(), x - 4 * sigma_), upper = std::upper_bound(x_.begin(), x_.end(), x + 4 * sigma_);
    for(std::vector<double>::const_iterator it = lower; it != upper; ++it)
    {
        const unsigned long i = it - x_.begin();
        check(*it == x_[i], "");

        const double k = kernel(x, (*it));
        res += k * y_[i];
        norm += k;
    }
    if(norm == 0)
        return 0;

    return res / norm;
}

double
GaussSmooth::kernel(double x, double x1) const
{
    check(sigmaSq_ > 0, "");
    const double diff = x - x1;
    const double diffSq = diff * diff;
    return std::exp(-diffSq / (2.0 * sigmaSq_));
}

class GaussSmooth2D : public Math::Function2<double, double, double>
{
public:
    inline GaussSmooth2D(const std::vector<double>& x1, const std::vector<double>& x2, const std::vector<std::vector<double> >& y, double sigma1, double sigma2 = 0);

    ~GaussSmooth2D() {}

    inline virtual double evaluate(double x1, double x2) const;

private:
    inline double kernel(double x1, double x2, double x1Prime, double x2Prime) const;

private:
    const std::vector<double> x1_, x2_;
    const std::vector<std::vector<double> > y_;

    const double sigmaSq1_, sigmaSq2_;
    const double sigma1_, sigma2_;
};

GaussSmooth2D::GaussSmooth2D(const std::vector<double>& x1, const std::vector<double>& x2, const std::vector<std::vector<double> >& y, double sigma1, double sigma2) : x1_(x1), x2_(x2), y_(y), sigmaSq1_(sigma1 * sigma1), sigmaSq2_((sigma2 == 0 ? sigma1 * sigma1 : sigma2 * sigma2)), sigma1_(sigma1), sigma2_(sigma2 == 0 ? sigma1 : sigma2)
{
    check(!x1_.empty(), "need to have at least 1 point");
    check(x1_.size() == y_.size(), "x1 must have the same size as y");
#ifdef CHECKS_ON
    for(unsigned long i = 0; i < y_.size(); ++i)
    {
        check(y_[i].size() == x2_.size(), "the elements of y must have the same size as x2");
    }
#endif

    check(sigma1 > 0, "invalid sigma1 = " << sigma1 << ", must be positive");
    check(sigma2 >= 0, "invalid sigma2 = " << sigma2 << ", must be positive or 0 to use sigma 1");
}

double
GaussSmooth2D::evaluate(double x1, double x2) const
{
    check(x1_.size() == y_.size(), "");

    double res = 0, norm = 0;
    const std::vector<double>::const_iterator lower1 = std::lower_bound(x1_.begin(), x1_.end(), x1 - 4 * sigma1_), upper1 = std::upper_bound(x1_.begin(), x1_.end(), x1 + 4 * sigma1_);
    for(std::vector<double>::const_iterator it1 = lower1; it1 != upper1; ++it1)
    {
        const int i = it1 - x1_.begin();
        check(x1_[i] == *it1, "");

        const std::vector<double>::const_iterator lower2 = std::lower_bound(x2_.begin(), x2_.end(), x2 - 4 * sigma2_), upper2 = std::upper_bound(x2_.begin(), x2_.end(), x2 + 4 * sigma2_);
        for(std::vector<double>::const_iterator it2 = lower2; it2 != upper2; ++it2)
        {
            const int j = it2 - x2_.begin();
            check(x2_[j] == *it2, "");
            const double k = kernel(x1, x2, *it1, *it2);
            res += k * y_[i][j];
            norm += k;
        }
    }
    if(norm > 0)
        return res / norm;

    return 0.0;
}

double
GaussSmooth2D::kernel(double x1, double x2, double x1Prime, double x2Prime) const
{
    check(sigmaSq1_ > 0 && sigmaSq2_ > 0, "");
    const double diff1 = x1 - x1Prime, diff2 = x2 - x2Prime;
    return std::exp(-diff1 * diff1 / (2.0 * sigmaSq1_) - diff2 * diff2 / (2.0 * sigmaSq2_));
}

} // namespace Math

#endif

