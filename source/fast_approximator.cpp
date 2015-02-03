#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <gaussian_process.hpp>
#include <fast_approximator.hpp>

#include <gmd.h>
#include <lavd.h>
#include <laslv.h>
#include <lavli.h>
#include <blas2pp.h>
#include <blas3pp.h>


FastApproximator::FastApproximator(int nPoints, int nData, unsigned long dataSize, const std::vector<std::vector<double> >& points, const std::vector<std::vector<double> >& data, int k) : knn_(NULL), k_(k), nPoints_(nPoints), nData_(nData), x_(k, nPoints + nPoints * nPoints + 1), xLin_(k, nPoints + 1), xT_(nPoints + nPoints * nPoints + 1, k), xTLin_(nPoints_ + 1, k), inv_(nPoints + nPoints * nPoints + 1, nPoints + nPoints * nPoints + 1), invLin_(nPoints_ + 1, nPoints_ + 1), prod_(nPoints + nPoints * nPoints + 1, k), prodLin_(nPoints_ + 1, k), pivot_(nPoints + nPoints * nPoints + 1), pivotLin_(nPoints_ + 1), pivot1_(nPoints), pointTransformed_(nPoints), sigma_(1), l_(1e-6), gaussK_(k, k), pivotK_(k)
{
    check(nPoints_ > 0, "");
    check(nData_ > 0, "");

    for(int i = 0; i < k; ++i)
    {
        x_(i, 0) = 1;
        xT_(0, i) = 1;

        xLin_(i, 0) = 1;
        xTLin_(0, i) = 1;
    }

    covariance_.resize(nPoints_, 2 * nPoints_ + 1);
    choleskyMat_.resize(nPoints_, nPoints_);
    v_.resize(nPoints_);
    w_.resize(nPoints_);

    indices_.resize(k_);
    dists_.resize(k_);

    reset(dataSize, points, data, true);
}

FastApproximator::~FastApproximator()
{
    check(knn_, "");
    delete knn_;
}

void
FastApproximator::reset(unsigned long dataSize, const std::vector<std::vector<double> >& points, const std::vector<std::vector<double> >& data, bool updateCovariance)
{
    output_screen("Fast Approximator learn with " << dataSize << " points." << std::endl);
    Timer timer("FAST APPROXIMATOR LEARN");
    timer.start();

    dataSize_ = dataSize;
    check(dataSize_ > 0, "");
    check(points.size() >= dataSize_, "");
    check(data.size() >= dataSize_, "");

    data_ = &data;

    if(updateCovariance)
    {
        std::vector<double> mean(nPoints_, 0);

        for(unsigned long i = 0; i < dataSize_; ++i)
        {
            check(points[i].size() >= nPoints_, "");

            for(int j = 0; j < nPoints_; ++j)
                mean[j] += points[i][j];
        }

        for(int i = 0; i < nPoints_; ++i)
            mean[i] /= double(dataSize_);

        for(int i = 0; i < nPoints_; ++i)
        {
            for(int j = 0; j < nPoints_; ++j)
            {
                double t1 = 0;
                for(unsigned long k = 0; k < dataSize_; ++k)
                    t1 += points[k][i] * points[k][j];

                covariance_(i, j) = (t1 - double(dataSize_) * mean[i] * mean[j]) / double(dataSize_ - 1);
            }
        }

        LaSymmBandMatFactorizeIP(covariance_);

        for(int i = 0; i < nPoints_; ++i)
        {
            for(int j = 0; j < nPoints_; ++j)
                choleskyMat_(i, j) = (j <= i ? covariance_(i, j) : 0.0);
        }

        LUFactorizeIP(choleskyMat_, pivot1_);
        LaLUInverseIP(choleskyMat_, pivot1_);
    }

    pointsTransformed_.resize(dataSize_);

    for(unsigned long i = 0; i < dataSize_; ++i)
    {
        pointsTransformed_[i].resize(nPoints_);

        for(int j = 0; j < nPoints_; ++j)
        {
            v_(j) = points[i][j];
        }

        Blas_Mat_Vec_Mult(choleskyMat_, v_, w_);
        for(int j = 0; j < nPoints_; ++j)
            pointsTransformed_[i][j] = w_(j);
    }

    if(!knn_)
        knn_ = new Math::KNearestNeighbors(nPoints_, pointsTransformed_, k_);
    else
        knn_->reset(pointsTransformed_);

    timer.end();
}

void
FastApproximator::findNearestNeighbors(const std::vector<double>& point, std::vector<double>* distances, std::vector<std::vector<double> >* nearestNeighbors, std::vector<unsigned long>* indices)
{
    check(point.size() == nPoints_, "");

    for(int i = 0; i < nPoints_; ++i)
    {
        v_(i) = point[i];
    }
    Blas_Mat_Vec_Mult(choleskyMat_, v_, w_);
    for(int i = 0; i < nPoints_; ++i)
        pointTransformed_[i] = w_(i);

    check(knn_, "");
    knn_->search(pointTransformed_, &indices_, &dists_);

    if(dists_[0] == 0)
    {
        output_screen("FOUND DISTANCE = 0" << std::endl);
        for(int i = 0; i < nPoints_; ++i)
        {
            //output_screen("\t" << pointTransformed_[i] << "\t" << pointsTransformed_[indices_[0]][i] << std::endl);
        }
    }

    if(distances)
    {
        distances->resize(k_);
        for(int i = 0; i < k_; ++i)
            (*distances)[i] = std::sqrt(dists_[i]);
    }

    if(nearestNeighbors)
    {
        nearestNeighbors->resize(k_);
        for(int i = 0; i < k_; ++i)
        {
            const unsigned long index = indices_[i];
            (*nearestNeighbors)[i] = pointsTransformed_[index];
            for(int j = 0; j < nPoints_; ++j)
                (*nearestNeighbors)[i][j] -= pointTransformed_[j];
        }
    }

    if(indices)
    {
        *indices = indices_;
    }
}

void
FastApproximator::getApproximationGaussianProcess(std::vector<double>& val, std::vector<double>& error)
{
    val.resize(nData_);
    error.resize(nData_);

    if(std::sqrt(dists_[0]) < 1e-7)
    {
        //output_screen("FOUND distance = " << dists_[0] << std::endl);
        const unsigned long index = indices_[0];
        for(int i = 0; i < nData_; ++i)
        {
            val[i] = (*data_)[index][i];
            error[i] = 0;
        }

        return;
    }

    std::vector<std::vector<double> > gaussPoints;
    for(int i = 0; i < k_; ++i)
        gaussPoints.push_back(pointsTransformed_[indices_[i]]);

    for(int i = 0; i < nData_; ++i)
    {
        std::vector<double> gaussData;
        for(int j = 0; j < k_; ++j)
            gaussData.push_back((*data_)[indices_[j]][i]);

        Math::GaussianProcess gp(nPoints_);
        gp.set(gaussPoints, gaussData, true);

        std::vector<double> mean;
        LaGenMatDouble covariance;

        std::vector<std::vector<double> > gaussInput(1, pointTransformed_);

        gp.calculate(gaussInput, mean, covariance);

        val[i] = mean[0];
        error[i] = std::sqrt(covariance(0, 0));
    }

    /*
    for(int i = 0; i < k_; ++i)
    {
        for(int j = i; j < k_; ++j)
        {
            const double c = cov(pointsTransformed_[indices_[i]], pointsTransformed_[indices_[j]]);
            gaussK_(i, j) = c;
            gaussK_(j, i) = c;

            //output_screen(i << ' ' << j << ' ' << c << std::endl);
        }
    }

    // hack
    for(int i = 0; i < k_; ++i)
        gaussK_(i, i) *= 1.00001;

    LUFactorizeIP(gaussK_, pivotK_);
    LaLUInverseIP(gaussK_, pivotK_);

    kStar_.resize(k_);
    for(int i = 0; i < k_; ++i)
        kStar_[i] = cov(dists_[i]);

    kKInv_.resize(k_);

    for(int i = 0; i < k_; ++i)
    {
        kKInv_[i] = 0;
        for(int j = 0; j < k_; ++j)
            kKInv_[i] += kStar_[j] * gaussK_(j, i);
    }

    double e = cov(0);
    for(int i = 0; i < k_; ++i)
        e -= kKInv_[i] * kStar_[i];

    check(e >= 0, "");
    e = std::sqrt(e);

    for(int i = 0; i < nData_; ++i)
    {
        error[i] = e;
        val[i] = 0;
        for(int j = 0; j < k_; ++j)
            val[i] += kKInv_[j] * (*data_)[indices_[j]][i];
    }
    */
}

void
FastApproximator::getApproximation(std::vector<double>& val, InterpolationMethod method)
{
    check(method >= 0 && method < INTERPOLATION_METHOD_MAX, "");

    val.resize(nData_);

    if(std::sqrt(dists_[0]) < 1e-7)
    {
        //output_screen("FOUND distance = " << dists_[0] << std::endl);
        const unsigned long index = indices_[0];
        for(int i = 0; i < nData_; ++i)
            val[i] = (*data_)[index][i];

        return;
    }

    //Timer t1("MATRIX OPERATIONS");
    //t1.start();

    std::vector<double> weights(k_);
    for(int i = 0; i < k_; ++i)
        weights[i] = 1.0 / std::sqrt(dists_[i]);
    
    for(int i = 0; i < k_; ++i)
    {
        const unsigned long index = indices_[i];
        for(int j = 0; j < nPoints_; ++j)
        {
            const double x = pointsTransformed_[index][j] - pointTransformed_[j];
            switch(method)
            {
            case LINEAR_INTERPOLATION:
                xLin_(i, j + 1) = x * weights[i];
                xTLin_(j + 1, i) = x;
                break;
            case QUADRATIC_INTERPOLATION:
                x_(i, j + 1) = x * weights[i];
                xT_(j + 1, i) = x;
                for(int l = 0; l < nPoints_; ++l)
                {
                    const double y = pointsTransformed_[index][l] - pointTransformed_[l];
                    x_(i, nPoints_ + 1 + j * nPoints_ + l) = x * y * weights[i];
                    xT_(nPoints_ + 1 + j * nPoints_ + l, i) = x * y;
                }
                break;
            default:
                check(false, "");
                break;
            }

        }

        switch(method)
        {
        case LINEAR_INTERPOLATION:
            xLin_(i, 0) = weights[i];
            break;
        case QUADRATIC_INTERPOLATION:
            x_(i, 0) = weights[i];
            break;
        default:
            check(false, "");
            break;
        }
    }

    switch(method)
    {
    case LINEAR_INTERPOLATION:
        Blas_Mat_Mat_Mult(xTLin_, xLin_, invLin_, false);

        for(int i = 0; i < invLin_.size(0); ++i)
            invLin_(i, i) += 1e-5;

        LUFactorizeIP(invLin_, pivotLin_);
        LaLUInverseIP(invLin_, pivotLin_);
        Blas_Mat_Mat_Mult(invLin_, xTLin_, prodLin_, false);
        break;

    case QUADRATIC_INTERPOLATION:
        Blas_Mat_Mat_Mult(xT_, x_, inv_, false);
        for(int i = 0; i < inv_.size(0); ++i)
            inv_(i, i) += 1e-5;

        LUFactorizeIP(inv_, pivot_);
        LaLUInverseIP(inv_, pivot_);
        Blas_Mat_Mat_Mult(inv_, xT_, prod_, false);
        break;

    default:
        check(false, "");
        break;
    }


    //t1.end();

    //Timer t2("PARAMETER ESTIMATION");
    //t2.start();

    for(int i = 0; i < nData_; ++i)
    {
        double res = 0;
        for(int j = 0; j < k_; ++j)
        {
            const unsigned long index = indices_[j];
            const double y = (*data_)[index][i] * weights[j];

            switch(method)
            {
            case LINEAR_INTERPOLATION:
                res += prodLin_(0, j) * y;
                break;
            case QUADRATIC_INTERPOLATION:
                res += prod_(0, j) * y;
                break;
            default:
                check(false, "");
                break;
            }
        }

        val[i] = res;
    }

    //t2.end();
}

void
FastApproximator::approximate(const std::vector<double>& point, std::vector<double>& val, InterpolationMethod method, std::vector<double>* distances, std::vector<std::vector<double> >* nearestNeighbors, std::vector<unsigned long>* indices)
{
    findNearestNeighbors(point, distances, nearestNeighbors, indices);
    getApproximation(val, method);
}

