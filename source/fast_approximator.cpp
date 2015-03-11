#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <matrix_impl.hpp>
#include <fast_approximator.hpp>

FastApproximator::FastApproximator(int nPoints, int nData, unsigned long dataSize, const std::vector<std::vector<double> >& points, const std::vector<std::vector<double> >& data, int k) : knn_(NULL), k_(k), nPoints_(nPoints), nData_(nData), x_(k, nPoints + nPoints * nPoints + 1), xLin_(k, nPoints + 1), xT_(nPoints + nPoints * nPoints + 1, k), xTLin_(nPoints_ + 1, k), inv_(nPoints + nPoints * nPoints + 1, nPoints + nPoints * nPoints + 1), invLin_(nPoints_ + 1, nPoints_ + 1), prod_(nPoints + nPoints * nPoints + 1, k), prodLin_(nPoints_ + 1, k), pointTransformed_(nPoints), sigma_(1), l_(1e-6)
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

    covariance_.resize(nPoints_);
    choleskyMat_.resize(nPoints_, nPoints_);
    v_.resize(nPoints_, 1);
    w_.resize(nPoints_, 1);

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
            for(int j = i; j < nPoints_; ++j)
            {
                double t1 = 0;
                for(unsigned long k = 0; k < dataSize_; ++k)
                    t1 += points[k][i] * points[k][j];

                covariance_(i, j) = (t1 - double(dataSize_) * mean[i] * mean[j]) / double(dataSize_ - 1);
            }
        }

        covariance_.choleskyFactorize();

        for(int i = 0; i < nPoints_; ++i)
        {
            for(int j = 0; j < nPoints_; ++j)
                choleskyMat_(i, j) = (j <= i ? covariance_(i, j) : 0.0);
        }

        choleskyMat_.invert();
    }

    pointsTransformed_.resize(dataSize_);

    for(unsigned long i = 0; i < dataSize_; ++i)
    {
        pointsTransformed_[i].resize(nPoints_);

        for(int j = 0; j < nPoints_; ++j)
        {
            v_(j, 0) = points[i][j];
        }

        Math::Matrix<double>::multiplyMatrices(choleskyMat_, v_, &w_);
        for(int j = 0; j < nPoints_; ++j)
            pointsTransformed_[i][j] = w_(j, 0);
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
        v_(i, 0) = point[i];
    }
    Math::Matrix<double>::multiplyMatrices(choleskyMat_, v_, &w_);
    for(int i = 0; i < nPoints_; ++i)
        pointTransformed_[i] = w_(i, 0);

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
        Math::Matrix<double>::multiplyMatrices(xTLin_, xLin_, &invLin_);

        for(int i = 0; i < invLin_.rows(); ++i)
            invLin_(i, i) += 1e-5;

        invLin_.invert();
        Math::Matrix<double>::multiplyMatrices(invLin_, xTLin_, &prodLin_);
        break;

    case QUADRATIC_INTERPOLATION:
        Math::Matrix<double>::multiplyMatrices(xT_, x_, &inv_);
        for(int i = 0; i < inv_.rows(); ++i)
            inv_(i, i) += 1e-5;

        inv_.invert();
        Math::Matrix<double>::multiplyMatrices(inv_, xT_, &prod_);
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

