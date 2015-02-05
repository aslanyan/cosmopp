#ifndef COSMO_PP_K_NEAREST_NEIGHBORS_HPP
#define COSMO_PP_K_NEAREST_NEIGHBORS_HPP

#include <vector>

#include <ANN/ANN.h>

namespace Math
{

class KNearestNeighbors
{
public:
    KNearestNeighbors(int dim, const std::vector<std::vector<double> >& points, int k, double epsilon = 0);
    ~KNearestNeighbors();

    void reset(const std::vector<std::vector<double> >& points);

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

