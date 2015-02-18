#ifndef COSMO_PP_K_NEAREST_NEIGHBORS_HPP
#define COSMO_PP_K_NEAREST_NEIGHBORS_HPP

#include <vector>

#include <ANN/ANN.h>

namespace Math
{

/// This class can be used to find k nearest neighbors of a given point in a given set of points. The Euclidean distance is used.
class KNearestNeighbors
{
public:
    /// Constructor.
    /// \param dim The dimensionality of space.
    /// \param points A vector of the points (each of which should be a dim dimenson vector) that form the set to be searched in.
    /// \param k The number of nearest neighbors to find.
    /// \param epsilon The error tolerance. The default value is 0 which means there is no error tolerance.
    KNearestNeighbors(int dim, const std::vector<std::vector<double> >& points, int k, double epsilon = 0);

    /// Desctructor.
    ~KNearestNeighbors();

    /// Reset the set of points to search in.
    /// \param points A vector of the points (each of which should be a dim dimenson vector) that form the set to be searched in.
    void reset(const std::vector<std::vector<double> >& points);

    /// Find the k nearest neighbors of a given point.
    /// \param point The point to find the nearest neighbors of.
    /// \param indices A pointer to a vector where the indices of the nearest neighbors will be returned. The indices will be the same as in the original input set.
    /// \param distanceSquares A pointer to a vector where the distances to the nearest neighbors squared will be returned. Can be NULL (by default) in which case the distances are not returned.
    void search(const std::vector<double>& point, std::vector<unsigned long>* indices, std::vector<double>* distanceSquares = NULL);

private:
    void deallocate();

private:
    struct SingleDestructorObject
    {
        SingleDestructorObject() {}
        ~SingleDestructorObject() { annClose(); }
    };

private:
    ANNpointArray dataPts_;
    ANNkd_tree* kdTree_;
    ANNpoint queryPt_;
    ANNidxArray nnIdx_;
    ANNdistArray dists_;

    const int dim_;
    const int k_;
    const double epsilon_;
};

} // namespace Math

#endif

