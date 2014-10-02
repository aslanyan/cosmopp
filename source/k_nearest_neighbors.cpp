#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <timer.hpp>
#include <k_nearest_neighbors.hpp>

namespace Math
{

KNearestNeighbors::KNearestNeighbors(int dim, const std::vector<std::vector<double> >& points, int k, double epsilon) : dim_(dim), k_(k), epsilon_(epsilon), nnIdx_(NULL), dists_(NULL), kdTree_(NULL), dataPts_(NULL), queryPt_(NULL)
{
    check(k_ > 0, "invalid k = " << k_ << ", must be positive");
    check(epsilon_ >= 0, "invalid epsilon = " << epsilon_ << ", cannot be negative");

    check(dim > 0, "invalid dim = " << dim << ", must be positive");

    reset(points);

    static SingleDestructorObject sdo;
}

KNearestNeighbors::~KNearestNeighbors()
{
    deallocate();
}

void
KNearestNeighbors::deallocate()
{
    if(nnIdx_)
    {
        delete [] nnIdx_;
        nnIdx_ = NULL;
    }

    if(dists_)
    {
        delete [] dists_;
        dists_ = NULL;
    }

    if(kdTree_)
    {
        delete kdTree_;
        kdTree_ = NULL;
    }

    if(dataPts_)
    {
        annDeallocPts(dataPts_);
        dataPts_ = NULL;
    }

    if(queryPt_)
    {
        annDeallocPt(queryPt_);
        queryPt_ = NULL;
    }
}

void
KNearestNeighbors::reset(const std::vector<std::vector<double> >& points)
{
    deallocate();

    check(!points.empty(), "need at least one point");

    Timer timer("ANN LEARN");
    timer.start();

    const unsigned long size = points.size();

    dataPts_ = annAllocPts(size, dim_);
    queryPt_ = annAllocPt(dim_);

    nnIdx_ = new ANNidx[k_];
    dists_ = new ANNdist[k_];

    for(unsigned long i = 0; i < size; ++i)
    {
        check(points[i].size() >= dim_, "point i " << " has only " << points[i].size() << " elements, but the dimension is " << dim_);

        for(int j = 0; j < dim_; ++j)
        {
            dataPts_[i][j] = points[i][j];
        }
    }

    kdTree_ = new ANNkd_tree(dataPts_, size, dim_);
    timer.end();
}

void
KNearestNeighbors::search(const std::vector<double>& point, std::vector<unsigned long>* indices, std::vector<double>* distanceSquares)
{
    check(point.size() >= dim_, "");

    for(int i = 0; i < dim_; ++i)
    {
        queryPt_[i] = point[i];
    }

    //Timer timer("ANN SEARCH");
    //timer.start();
    kdTree_->annkSearch(queryPt_, k_, nnIdx_, dists_, epsilon_);
    //timer.end();
    
    indices->resize(k_);
    for(int i = 0; i < k_; ++i)
        (*indices)[i] = nnIdx_[i];

    if(distanceSquares)
    {
        distanceSquares->resize(k_);

        for(int i = 0; i < k_; ++i)
            (*distanceSquares)[i] = dists_[i];
    }
}

} // namespace Math

