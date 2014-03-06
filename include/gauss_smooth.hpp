#ifndef COSMO_PP_GAUSS_SMOOTH
#define COSMO_PP_GAUSS_SMOOTH

#include <vector>
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
};

GaussSmooth::GaussSmooth(const std::vector<double>& x, const std::vector<double>& y, double sigma) : x_(x), y_(y), sigmaSq_(sigma * sigma)
{
    check(!x_.empty(), "need to have at least 1 point");
    check(x_.size() == y_.size(), "the sizes of vectors x and y must match");
    check(sigma > 0, "invalid sigma = " << sigma << ", must be positive");
}

double
GaussSmooth::evaluate(double x) const
{
    check(x_.size() == y_.size(), "");
    double res = 0, norm = 0;
    for(int i = 0; i < x_.size(); ++i)
    {
        const double k = kernel(x, x_[i]);
        res += k * y_[i];
        norm += k;
    }
    check(norm > 0, "");
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
    inline GaussSmooth2D(const std::vector<double>& x1, const std::vector<double>& x2, const std::vector<double>& y, double sigma1, double sigma2 = 0);

    ~GaussSmooth2D() {}

    inline virtual double evaluate(double x1, double x2) const;

private:
    inline double kernel(double x1, double x2, double x1Prime, double x2Prime) const;

private:
    const std::vector<double> x1_, x2_, y_;
    const double sigmaSq1_, sigmaSq2_;
};

GaussSmooth2D::GaussSmooth2D(const std::vector<double>& x1, const std::vector<double>& x2, const std::vector<double>& y, double sigma1, double sigma2) : x1_(x1), x2_(x2), y_(y), sigmaSq1_(sigma1 * sigma1), sigmaSq2_((sigma2 == 0 ? sigma1 * sigma1 : sigma2 * sigma2))
{
    check(!x1_.empty(), "need to have at least 1 point");
    check(x1_.size() == x2_.size(), "the sizes of vectors x1 and x2 must match");
    check(x1_.size() == y_.size(), "the size of vector y must match the sizes of x1 and x2");
    check(sigma1 > 0, "invalid sigma1 = " << sigma1 << ", must be positive");
    check(sigma2 >= 0, "invalid sigma2 = " << sigma2 << ", must be positive or 0 to use sigma 1");
}

double
GaussSmooth2D::evaluate(double x1, double x2) const
{
    check(x1_.size() == x2_.size(), "");
    check(x1_.size() == y_.size(), "");

    double res = 0, norm = 0;
    for(int i = 0; i < x1_.size(); ++i)
    {
        const double k = kernel(x1, x2, x1_[i], x2_[i]);
        res += k * y_[i];
        norm += k;
    }
    check(norm > 0, "");
    return res / norm;
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

