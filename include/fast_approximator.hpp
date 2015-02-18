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

/// A class that can be used to find approximate values of a function at a given point by using a training set with exact input and output values of the function.
/// There are training input and output points of certain dimensions. We call the inputs of the function "points" and the outputs "data".
/// These training points are then used to find the output of a new input point by first finding k nearest neighbors of the point in the training set, then performing an interpolation between them.
/// The input space is linearly transformed to decorrelate the input parameters, i.e. the covariance matrix of the input parameters is calculated on the training set, the Cholesky decomposed, and this matrix is used to linearly transform the space to another basis where the input parameters are uncorrelated for the training set.
/// This is important because different input parameters may have different magnitudes, and there may be significant correlations between them. Since the algorithm uses k nearest neighbors for the approximation, it is important that different input parameters contribute to the distance equally.
/// This class is in fact a machine learning regression class.
class FastApproximator
{
public:
    /// The interpolation method.
    enum InterpolationMethod { LINEAR_INTERPOLATION = 0, QUADRATIC_INTERPOLATION, INTERPOLATION_METHOD_MAX };

public:
    /// Constructor.
    /// \param nPoints The dimensionality of the points space, i.e. the number of the input parameters.
    /// \param nData The dimensionality of the data space, i.e. the number of the output parameters.
    /// \param dataSize The number of data points to use.
    /// \param points A vector containing all of the input points. Each point should be a vector of dimension nPoints. There needs to be at least dataSize points here. If the size of this vector is larger than dataSize then only the first dataSize points will be used.
    /// \param data A vector containing all of the output points. Each point should be a vector of dimension nData. There needs to be at least dataSize points here. If the size of this vector is larger than dataSize then only the first dataSize points will be used. The indices of data should exactly match the indices of points.
    /// \param k The number of nearest neighbors to use in the approximation.
    FastApproximator(int nPoints, int nData, unsigned long dataSize, const std::vector<std::vector<double> >& points, const std::vector<std::vector<double> >& data, int k);

    /// Destructor.
    ~FastApproximator();

    /// Reset the training set.
    /// \param dataSize The number of data points to use.
    /// \param points A vector containing all of the input points. Each point should be a vector of dimension nPoints. There needs to be at least dataSize points here. If the size of this vector is larger than dataSize then only the first dataSize points will be used.
    /// \param data A vector containing all of the output points. Each point should be a vector of dimension nData. There needs to be at least dataSize points here. If the size of this vector is larger than dataSize then only the first dataSize points will be used. The indices of data should exactly match the indices of points.
    /// \param updateCovariance If this is set to true (by default) then the covariance matrix of the input parameters is recalculated for the new training set, and the linear transformation matrix is updated. It is important to keep in mind that if this step is performed then the distances to previously existing points will change. For example, if the training set is updated by just adding some new points and we want to keep the distances to the old points unchanged then this parameter should be set to false.
    void reset(unsigned long dataSize, const std::vector<std::vector<double> >& points, const std::vector<std::vector<double> >& data, bool updateCovariance = true);

    /// Find the nearest neighbors to a given point. This step needs to be always performed before calling getApproximation.
    /// \param point The input point.
    /// \param distances The distances to the nearest neighbors squared will be returned in this vector. This can be set to NULL if the distances are not needed (by default). Keep in mind that the distances are the Euclidean distances in a linearly transformed space where the input training parameters are decorrelated.
    /// \param nearestNeighbors The nearest neighbors RELATIVE to the input point in the linearly transformed space will be returned here. This can be set to NULL if not needed (by default).
    /// \param indices The indices of the nearest neighbors will be returned here. Set to NULL if not needed (by default).
    void findNearestNeighbors(const std::vector<double>& point, std::vector<double>* distances = NULL, std::vector<std::vector<double> >* nearestNeighbors = NULL, std::vector<unsigned long>* indices = NULL);

    /// Get the approximation of the output for the input point given to findNearestNeighbors. This function should be called after findNearestNeighbors.
    /// \param val The output will be returned here.
    /// \param method The interpolation method to be used.
    void getApproximation(std::vector<double>& val, InterpolationMethod method = QUADRATIC_INTERPOLATION);
    
    /// Experimental function.
    void getApproximationGaussianProcess(std::vector<double>& val, std::vector<double>& error);

    /// Find the approximate output for a given input point. This function is equivalent to calling findNearestNeighbors followed by getApproximation.
    /// \param point The input point.
    /// \param val The output will be returned here.
    /// \param method The interpolation method to be used.
    /// \param distances The distances to the nearest neighbors squared will be returned in this vector. This can be set to NULL if the distances are not needed (by default). Keep in mind that the distances are the Euclidean distances in a linearly transformed space where the input training parameters are decorrelated.
    /// \param nearestNeighbors The nearest neighbors RELATIVE to the input point in the linearly transformed space will be returned here. This can be set to NULL if not needed (by default).
    /// \param indices The indices of the nearest neighbors will be returned here. Set to NULL if not needed (by default).
    void approximate(const std::vector<double>& point, std::vector<double>& val, InterpolationMethod = QUADRATIC_INTERPOLATION, std::vector<double>* distances = NULL, std::vector<std::vector<double> >* nearestNeighbors = NULL, std::vector<unsigned long>* indices = NULL);

    /// Get the dimensionality of the input space.
    int nPoints() const { return nPoints_; }

    /// Get the dimensionality of the output space.
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

