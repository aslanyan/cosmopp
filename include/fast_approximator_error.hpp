#ifndef COSMO_PP_FAST_APPROXIMATOR_ERROR_HPP
#define COSMO_PP_FAST_APPROXIMATOR_ERROR_HPP

#include <fast_approximator.hpp>
#include <markov_chain.hpp>

template<typename Function>
class FastApproximatorError
{
public:
    enum ErrorMethod { MIN_DISTANCE = 0, AVG_DISTANCE, AVG_INV_DISTANCE, SUM_DISTANCE, LIN_QUAD_DIFF, ERROR_METHOD_MAX };

public:
    FastApproximatorError(FastApproximator& fa, const std::vector<std::vector<double> >& testPoints, const std::vector<std::vector<double> >& testData, unsigned long begin, unsigned long end, Function& f, ErrorMethod method = AVG_DISTANCE, double precision = 1.0);
    ~FastApproximatorError();

    void reset(const std::vector<std::vector<double> >& testPoints, const std::vector<std::vector<double> >& testData, unsigned long begin, unsigned long end);

    bool approximate(const std::vector<double>& point, std::vector<double>& val);

    void setPrecision(double p) { check(p > 0, "invalid precision " << p); precision_ = p; }

    Posterior1D* getPosterior() { return posterior_; }

private:
    double evaluateError();

private:
    Function& f_;
    FastApproximator& fa_;
    ErrorMethod method_;

    Posterior1D* posterior_;
    bool posteriorGood_;

    std::vector<double> val_;
    std::vector<double> linVal_;

    std::vector<double>* distances_;
    std::vector<std::vector<double> >* nearestNeighbors_;

    std::vector<double> distanceSum_;

    double precision_;
};

class BasicFAErrorFunctionAvg
{
public:
    double evaluate(const std::vector<double>& v) const
    {
        return v[0];

        double res = 0;
        for(int i = 0; i < v.size(); ++i)
            res += v[i];

        return res / v.size();
    }
};



template<typename Function>
FastApproximatorError<Function>::FastApproximatorError(FastApproximator& fa, const std::vector<std::vector<double> >& testPoints, const std::vector<std::vector<double> >& testData, unsigned long begin, unsigned long end, Function& f, ErrorMethod method, double precision) : fa_(fa), method_(method), posterior_(NULL), distances_(NULL), nearestNeighbors_(NULL), val_(fa.nData()), linVal_(fa.nData()), f_(f), precision_(precision), posteriorGood_(false)
{
    check(precision_ > 0, "invalid precision " << precision_);

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

template<typename Function>
FastApproximatorError<Function>::~FastApproximatorError()
{
    check(posterior_, "");
    delete posterior_;

    if(distances_)
        delete distances_;

    if(nearestNeighbors_)
        delete nearestNeighbors_;
}

template<typename Function>
void
FastApproximatorError<Function>::reset(const std::vector<std::vector<double> >& testPoints, const std::vector<std::vector<double> >& testData, unsigned long begin, unsigned long end)
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

    ProgressMeter meter(end - begin);
    unsigned long goodCount = 0;
    for(unsigned long i = begin; i < end; ++i)
    {
        fa_.approximate(testPoints[i], val_, FastApproximator::QUADRATIC_INTERPOLATION, distances_, nearestNeighbors_);

        if(method_ == LIN_QUAD_DIFF)
            fa_.getApproximation(linVal_, FastApproximator::LINEAR_INTERPOLATION);

        const double estimatedError = evaluateError();
        const double correctError = std::abs(f_.evaluate(testData[i]) - f_.evaluate(val_));

        //output_screen1("Estimated error = " << estimatedError << std::endl << "Correct error = " << correctError << std::endl << testData[i][100] << '\t' << val_[100] << std::endl);
        
        if(estimatedError == 0)
        {
            check(correctError == 0, "");
        }
        else
        {
            const double ratio = correctError / estimatedError;
            posterior_->addPoint(ratio, 1, 1);
            ++goodCount;
        }
        meter.advance();
    }

    if(goodCount >= 100)
    {
        posterior_->generate();
        posteriorGood_ = true;
        output_screen1("Posterior 1 sigma is: " << posterior_->get1SigmaUpper() << std::endl);
    }
    else
        posteriorGood_ = false;
    posterior_->writeIntoFile("fast_approximator_error_ratio.txt");

    t.end();
}

template<typename Function>
double
FastApproximatorError<Function>::evaluateError()
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

template<typename Function>
bool
FastApproximatorError<Function>::approximate(const std::vector<double>& point, std::vector<double>& val)
{
    fa_.findNearestNeighbors(point, distances_, nearestNeighbors_);
    if(method_ == LIN_QUAD_DIFF)
    {
        fa_.getApproximation(val_);
        fa_.getApproximation(linVal_, FastApproximator::LINEAR_INTERPOLATION);
    }
    const double e = evaluateError();
    if(!posteriorGood_)
    {
        if(e == 0)
        {
            output_screen1("Error = " << 0 << std::endl);
            fa_.getApproximation(val);
            return true;
        }
        return false;
    }

    const double estimatedError = e * posterior_->get1SigmaUpper(); //generateSample();

    output_screen1("Error = " << estimatedError << std::endl);

    if(estimatedError > precision_)
        return false;

    if(method_ == LIN_QUAD_DIFF)
        val = val_;
    else
        fa_.getApproximation(val);
    return true;
}
#endif

