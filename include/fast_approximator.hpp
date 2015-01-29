#ifndef COSMO_PP_FAST_APPROXIMATOR_HPP
#define COSMO_PP_FAST_APPROXIMATOR_HPP

#include <vector>

#include <macros.hpp>
#include <k_nearest_neighbors.hpp>
#include <timer.hpp>
#include <progress_meter.hpp>

#include <lavd.h>
#include <gmd.h>
#include <lavli.h>
#include <sybmd.h>
#include <sybfd.h>

#include <ANN/ANN.h>

class FastApproximator
{
public:
    enum InterpolationMethod { LINEAR_INTERPOLATION = 0, QUADRATIC_INTERPOLATION, INTERPOLATION_METHOD_MAX };
public:
    FastApproximator(int nPoints, int nData, unsigned long dataSize, const std::vector<std::vector<double> >& points, const std::vector<std::vector<double> >& data, int k);
    ~FastApproximator();

    void reset(unsigned long dataSize, const std::vector<std::vector<double> >& points, const std::vector<std::vector<double> >& data, bool updateCovariance = true);

    void findNearestNeighbors(const std::vector<double>& point, std::vector<double>* distances = NULL, std::vector<std::vector<double> >* nearestNeighbors = NULL, std::vector<unsigned long>* indices = NULL);
    void getApproximation(std::vector<double>& val, InterpolationMethod method = QUADRATIC_INTERPOLATION);
    void getApproximationGaussianProcess(std::vector<double>& val, std::vector<double>& error);

    void approximate(const std::vector<double>& point, std::vector<double>& val, InterpolationMethod = QUADRATIC_INTERPOLATION, std::vector<double>* distances = NULL, std::vector<std::vector<double> >* nearestNeighbors = NULL, std::vector<unsigned long>* indices = NULL);

    int nPoints() const { return nPoints_; }
    int nData() const { return nData_; }

private:
    inline double cov(double d) const { return sigma_ * std::exp(-d / (2 * l_)); }
    inline double cov(const std::vector<double>& x, const std::vector<double>& y) const
    {
        double d = 0;
        check(x.size() == nPoints_, "");
        check(y.size() == nPoints_, "");

        for(int i = 0; i < nPoints_; ++i)
        {
            const double delta = x[i] - y[i];
            d += delta * delta;
        }

        return cov(d);
    }

private:
    const int k_;
    Math::KNearestNeighbors* knn_;

    const std::vector<std::vector<double > >* data_;

    std::vector<std::vector<double> > pointsTransformed_;

    std::vector<double> pointTransformed_;

    unsigned long dataSize_;
    int nPoints_;
    int nData_;

    std::vector<unsigned long> indices_;
    std::vector<double> dists_;

    double sigma_, l_;

    LaSymmBandMatDouble covariance_;
    LaGenMatDouble choleskyMat_;
    LaVectorLongInt pivot1_;
    LaVectorDouble v_, w_;

    LaGenMatDouble x_, xLin_;
    LaGenMatDouble xT_, xTLin_;
    LaGenMatDouble inv_, invLin_;
    LaGenMatDouble prod_, prodLin_;
    LaVectorLongInt pivot_, pivotLin_;

    LaGenMatDouble gaussK_;
    LaVectorLongInt pivotK_;
    std::vector<double> kStar_;
    std::vector<double> kKInv_;
};

#endif

