#ifndef COSMO_PP_FAST_APPROXIMATOR_ERROR_HPP
#define COSMO_PP_FAST_APPROXIMATOR_ERROR_HPP

#include <fstream>

#include <function.hpp>
#include <fast_approximator.hpp>
#include <markov_chain.hpp>

/// A class for error evaluation for the FastApproximator class.
class FastApproximatorError
{
public:
    /// Error evaluation method.
    enum ErrorMethod { MIN_DISTANCE = 0, AVG_DISTANCE, AVG_INV_DISTANCE, SUM_DISTANCE, LIN_QUAD_DIFF, ERROR_METHOD_MAX };

    /// What to use in the posteriors to make a decision whether or not the error is acceptable. For example, for ONE_SIGMA, the one sigma (68%) upper bound of the error probability distribution will be compared to the precision to decide if the error is small enough to be acceptable.
    enum DecisionMethod { ONE_SIGMA = 0, TWO_SIGMA, SQRT_VAR, DECISION_METHOD_MAX };

public:
    /// Constructor.
    /// \param fa A reference to the fast approximator being used.
    /// \param testPoints The input points of the test set (the cross-validation set) that are used to evaluate the error distribution.
    /// \param testValues The output points of the test set (the cross-validation set) that are used to evaluate the error distribution. Must have the same size as testPoints.
    /// \param begin The starting index of the testPoints and testValues vectors. To use all of the points set begin to 0.
    /// \param end The index after the last point to be used. To use all of the points set end to the size of testPoints.
    /// \param f A function used to evaluate the error for. This function takes as an input the output of the fast approximator and gives as an output a single value for which the error will be evaluated.
    /// \param method The method to be used to evaluate the error.
    /// \param precision This is the error threshold. The error value will be acceptable if it's smaller than precision.
    /// \param dm The decision method, i.e. what property of the error probability distribution to use to compare to precision.
    FastApproximatorError(FastApproximator& fa, const std::vector<std::vector<double> >& testPoints, const std::vector<std::vector<double> >& testValues, unsigned long begin, unsigned long end, const Math::RealFunctionMultiDim& f, ErrorMethod method = AVG_DISTANCE, double precision = 1.0, DecisionMethod dm = TWO_SIGMA);

    /// Destructor.
    ~FastApproximatorError();

    /// Reset the test set (the cross-validation set).
    /// \param testPoints The input points of the test set (the cross-validation set) that are used to evaluate the error distribution.
    /// \param testValues The output points of the test set (the cross-validation set) that are used to evaluate the error distribution. Must have the same size as testPoints.
    /// \param begin The starting index of the testPoints and testValues vectors. To use all of the points set begin to 0.
    /// \param end The index after the last point to be used. To use all of the points set end to the size of testPoints.
    void reset(const std::vector<std::vector<double> >& testPoints, const std::vector<std::vector<double> >& testValues, unsigned long begin, unsigned long end);

    /// Approximate function.
    /// \param point The input point at which the approximation needs to be done.
    /// \param val The approximated result is returned here.
    /// \param error1Sigma If specified (i.e. not NULL), the one sigma upper bound of the absolute error probability distribution will be returned here.
    /// \param error2Sigma If specified (i.e. not NULL), the two sigma upper bound of the absolute error probability distribution will be returned here.
    /// \param errorMean If specified (i.e. not NULL), the mean of the error probability distribution will be returned here.
    /// \param errorVar If specified (i.e. not NULL), the variance of the error probability distribution will be returned here.
    bool approximate(const std::vector<double>& point, std::vector<double>& val, double *error1Sigma = NULL, double *error2Sigma = NULL, double *errorMean = NULL, double *errorVar = NULL);

    /// Set the precision.
    /// \param precision This is the error threshold. The error value will be acceptable if it's smaller than precision.
    /// \param dm The decision method, i.e. what property of the error probability distribution to use to compare to precision.
    void setPrecision(double p, DecisionMethod dm = TWO_SIGMA) { check(p > 0, "invalid precision " << p); check(dm >= 0 && dm < DECISION_METHOD_MAX, ""); precision_ = p; }

    /// Get the error probability distribution.
    Posterior1D* getDistrib() { return posterior_; }

private:
    double evaluateError();

private:
    const Math::RealFunctionMultiDim& f_;
    FastApproximator& fa_;
    ErrorMethod method_;

    Posterior1D* posterior_;
    bool posteriorGood_;
    double mean_, var_;

    std::vector<double> val_;
    std::vector<double> linVal_;

    std::vector<double>* distances_;
    std::vector<std::vector<double> >* nearestNeighbors_;

    std::vector<double> distanceSum_;

    std::vector<double> ge_;

    double precision_;
    DecisionMethod decMethod_;
};

class BasicFAErrorFunctionAvg : public Math::RealFunctionMultiDim
{
public:
    double evaluate(const std::vector<double>& v) const
    {
        double res = 0;
        for(int i = 0; i < v.size(); ++i)
            res += v[i];

        return res / v.size();
    }
};

#endif

