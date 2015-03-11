#include <sstream>
#include <string>
#include <iomanip>

#include <exception_handler.hpp>
#include <fast_approximator_error.hpp>

FastApproximatorError::FastApproximatorError(FastApproximator& fa, const std::vector<std::vector<double> >& testPoints, const std::vector<std::vector<double> >& testData, unsigned long begin, unsigned long end, const Math::RealFunctionMultiDim& f, ErrorMethod method, double precision, DecisionMethod dm) : fa_(fa), method_(method), posterior_(NULL), distances_(NULL), nearestNeighbors_(NULL), val_(fa.nData()), linVal_(fa.nData()), f_(f), precision_(precision), decMethod_(dm), posteriorGood_(false), mean_(0), var_(0)
{
    check(precision_ > 0, "invalid precision " << precision_);
    check(decMethod_ >= 0 && decMethod_ < DECISION_METHOD_MAX, "");

    check(method_ >= 0 && method_ < ERROR_METHOD_MAX, "invalid method");

    switch(method_)
    {
    case MIN_DISTANCE:
    case AVG_DISTANCE:
    case AVG_INV_DISTANCE:
        distances_ = new std::vector<double>;
        break;
    case SUM_DISTANCE:
        nearestNeighbors_ = new std::vector<std::vector<double> >;
        break;
    case LIN_QUAD_DIFF:
        break;
    default:
        check(false, "");
        break;
    }

    reset(testPoints, testData, begin, end);
}

FastApproximatorError::~FastApproximatorError()
{
    check(posterior_, "");
    delete posterior_;

    if(distances_)
        delete distances_;

    if(nearestNeighbors_)
        delete nearestNeighbors_;
}

void
FastApproximatorError::reset(const std::vector<std::vector<double> >& testPoints, const std::vector<std::vector<double> >& testData, unsigned long begin, unsigned long end)
{
    check(testPoints.size() >= end, "");
    check(testData.size() >= end, "");

    if(end == begin)
    {
        posteriorGood_ = false;
        return;
    }

    check(end > begin, "");

    Timer t("ERROR EVALUATION");
    t.start();

    if(posterior_)
        delete posterior_;

    posterior_ = new Posterior1D;
    mean_ = 0;
    var_ = 0;
    double mean2 = 0;
    unsigned long goodCount = 0;

    ProgressMeter meter(end - begin);
    for(unsigned long i = begin; i < end; ++i)
    {
        fa_.findNearestNeighbors(testPoints[i], distances_, nearestNeighbors_);
        fa_.getApproximation(val_, FastApproximator::QUADRATIC_INTERPOLATION);

        if(method_ == LIN_QUAD_DIFF)
            fa_.getApproximation(linVal_, FastApproximator::LINEAR_INTERPOLATION);

        const double estimatedError = evaluateError();
        const double correctError = f_.evaluate(testData[i]) - f_.evaluate(val_);

        if(estimatedError == 0)
        {
            check(correctError == 0, "");
        }
        else
        {
            const double ratio = correctError / estimatedError;
            posterior_->addPoint(std::abs(ratio), 1, 1);
            mean_ += ratio;
            mean2 += ratio * ratio;
            ++goodCount;
        }
        meter.advance();
    }

    if(goodCount >= 100)
    {
        posterior_->generate();
        posteriorGood_ = true;
        output_screen1("Posterior 1 sigma is: " << posterior_->get1SigmaUpper() << std::endl);
        output_screen1("Posterior 2 sigma is: " << posterior_->get2SigmaUpper() << std::endl);
        posterior_->writeIntoFile("fast_approximator_error_ratio.txt");

        mean_ /= goodCount;
        mean2 /= goodCount;

        var_ = mean2 - mean_ * mean_;
        check(var_ >= 0, "");

        output_screen("Error ratio distrib = " << mean_ << " +/- " << std::sqrt(var_) << std::endl);
    }
    else
    {
        posteriorGood_ = false;
        mean_ = 0;
        var_ = 0;
    }

    t.end();
}

double
FastApproximatorError::evaluateError()
{
    check(method_ >= 0 && method_ < ERROR_METHOD_MAX, "invalid method");

    double avgDist = 0;

    double sumDist = 0;

    double y, linY;

    switch(method_)
    {
    case MIN_DISTANCE:
        check(distances_, "");
        check(!distances_->empty(), "");
        return (*distances_)[0];

    case AVG_DISTANCE:
        check(distances_, "");
        check(!distances_->empty(), "");

        for(int i = 0; i < distances_->size(); ++i)
            avgDist += (*distances_)[i];

        avgDist /= distances_->size();
        return avgDist;

    case AVG_INV_DISTANCE:
        check(distances_, "");
        check(!distances_->empty(), "");

        for(int i = 0; i < distances_->size(); ++i)
        {
            if((*distances_)[i] == 0)
                return 0;

            avgDist += 1.0 / (*distances_)[i];
        }

        avgDist = distances_->size() / avgDist;

        return avgDist;

    case SUM_DISTANCE:
        check(nearestNeighbors_, "");
        check(!nearestNeighbors_->empty(), "");

        distanceSum_.resize((*nearestNeighbors_)[0].size(), 0);

        for(int i = 0; i < distanceSum_.size(); ++i)
            distanceSum_[i] = 0;

        for(int i = 0; i < nearestNeighbors_->size(); ++i)
        {
            check((*nearestNeighbors_)[i].size() == distanceSum_.size(), "");
            for(int j = 0; j < distanceSum_.size(); ++j)
                distanceSum_[j] += (*nearestNeighbors_)[i][j];
        }

        for(int i = 0; i < distanceSum_.size(); ++i)
            sumDist += distanceSum_[i] * distanceSum_[i];

        sumDist = std::sqrt(sumDist);

        return sumDist;

    case LIN_QUAD_DIFF:
        y = f_.evaluate(val_);
        linY = f_.evaluate(linVal_);
        //output_screen("Value = " << y << ", linear value = " << linY << std::endl);
        return std::abs(y - linY);
        
    default:
        check(false, "");
        return 0;
    }
}

bool
FastApproximatorError::approximate(const std::vector<double>& point, std::vector<double>& val, double *error1Sigma, double *error2Sigma, double *errorMean, double *errorVar)
{
    fa_.findNearestNeighbors(point, distances_, nearestNeighbors_);
    if(method_ == LIN_QUAD_DIFF)
    {
        fa_.getApproximation(val_);
        fa_.getApproximation(linVal_, FastApproximator::LINEAR_INTERPOLATION);
    }
    const double e = evaluateError();

    double estimatedError1 = 1e10, estimatedError2 = 1e10, estMean = 0, estVar = 1e20;
    if(posteriorGood_)
    {
        estimatedError1 = e * posterior_->get1SigmaUpper();
        estimatedError2 = e * posterior_->get2SigmaUpper();
        estMean = e * mean_;
        estVar = e * e * var_;
    }

    if(error1Sigma) *error1Sigma = estimatedError1;
    if(error2Sigma) *error2Sigma = estimatedError2;
    if(errorMean) *errorMean = estMean;
    if(errorVar) *errorVar = estVar;

    double errorChoice;
    switch(decMethod_)
    {
    case ONE_SIGMA:
        errorChoice = estimatedError1;
        break;
    case TWO_SIGMA:
        errorChoice = estimatedError2;
        break;
    case SQRT_VAR:
        errorChoice = std::sqrt(estVar);
        break;
    default:
        check(false, "");
        break;
    }

    if(errorChoice > precision_)
        return false;

    if(method_ == LIN_QUAD_DIFF)
        val = val_;
    else
        fa_.getApproximation(val);
    return true;
}
