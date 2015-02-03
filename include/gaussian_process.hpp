#ifndef COSMO_PP_GAUSSIAN_PROCESS_HPP
#define COSMO_PP_GAUSSIAN_PROCESS_HPP

#include <vector>

#include <lavd.h>
#include <gmd.h>
#include <lavli.h>
#include <sybmd.h>
#include <sybfd.h>

namespace Math
{

class GaussianProcess
{
public:
    enum CovarianceType { SQUARED_EXP = 0, COVARIANCE_TYPE_MAX };

public:
    GaussianProcess(int dim, CovarianceType covType = SQUARED_EXP);

    void set(const std::vector<std::vector<double> >& x, const std::vector<double>& y, bool calculateParamsByMaxLike = false);

    void calculate(const std::vector<std::vector<double> >& input, std::vector<double>& mean, LaGenMatDouble& covariance) const;

    void setParams(const std::vector<double>& params);

private:
    double cov(const std::vector<double>& x1, const std::vector<double>& x2) const;

private:
    int dim_;
    CovarianceType covType_;

    LaGenMatDouble k_;
    std::vector<double> y_;
    std::vector<std::vector<double> > x_;

    // squared exponential covariance params
    double sigma2_;
    std::vector<double> l2_;
};

} // namespace Math

#endif

