#ifndef COSMO_PP_FAST_APPROXIMATOR_ERROR_HPP
#define COSMO_PP_FAST_APPROXIMATOR_ERROR_HPP

#include <fstream>

#include <function.hpp>
#include <fast_approximator.hpp>
#include <markov_chain.hpp>

class FastApproximatorError
{
public:
    enum ErrorMethod { MIN_DISTANCE = 0, GAUSS_PROCESS, AVG_DISTANCE, AVG_INV_DISTANCE, SUM_DISTANCE, LIN_QUAD_DIFF, ERROR_METHOD_MAX };

public:
    FastApproximatorError(FastApproximator& fa, const std::vector<std::vector<double> >& testPoints, const std::vector<std::vector<double> >& testData, unsigned long begin, unsigned long end, const Math::RealFunctionMultiDim& f, ErrorMethod method = AVG_DISTANCE, double precision = 1.0);
    ~FastApproximatorError();

    void reset(const std::vector<std::vector<double> >& testPoints, const std::vector<std::vector<double> >& testData, unsigned long begin, unsigned long end);

    bool approximate(const std::vector<double>& point, std::vector<double>& val);

    void setPrecision(double p) { check(p > 0, "invalid precision " << p); precision_ = p; }

    Posterior1D* getPosterior() { return posterior_; }

    void logIntoFile(const char* fileNameBase);

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

    double gaussError_;
    std::vector<double> ge_;

    double precision_;

    std::ofstream outLog_;
    bool logFileOpen_;
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

