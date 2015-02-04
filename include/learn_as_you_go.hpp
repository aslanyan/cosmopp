#ifndef COSMO_PP_LEARN_AS_YOU_GO_HPP
#define COSMO_PP_LEARN_AS_YOU_GO_HPP

#include <fstream>
#include <vector>
#include <string>
#include <map>

#include <macros.hpp>
#include <function.hpp>
#include <random.hpp>
#include <fast_approximator.hpp>
#include <fast_approximator_error.hpp>

class LearnAsYouGo
{
public:
    LearnAsYouGo(int nPoints, int nData, const Math::RealFunctionMultiToMulti& f, const Math::RealFunctionMultiDim& errorFunc, unsigned long minCount = 1000, double precision = 0.1, const char* fileName = "");

    ~LearnAsYouGo();
    
    void evaluate(const std::vector<double>& x, std::vector<double>* res);

    void setPrecision(double p);

    void writeIntoFile(const char* fileName) const;
    bool readFromFile(const char* fileName);

    void logIntoFile(const char* fileNameBase);

    unsigned long getTotalCount() const { return totalCount_; }
    unsigned long getSuccessfulCount() const { return successfulCount_; }

private:
    void construct();
    void randomizeErrorSet();

    void constructFast();

    void log();

    void actual(const std::vector<double>& x, std::vector<double>* res);

    void communicate();
    void receive();

    void addDataPoint(const std::vector<double>& p, const std::vector<double>& d);

    void resetPointMap();

private:
    struct PointComp
    {
        bool operator()(const std::vector<double>& a, const std::vector<double>& b) const
        {
            check(a.size() == b.size(), "");
            check(!a.empty(), "");

            for(int i = 0; i < a.size(); ++i)
            {
                if(a[i] < b[i])
                    return true;
                if(a[i] > b[i])
                    return false;
            }

            return false;
        }
    };

private:
    int nPoints_, nData_;
    const Math::RealFunctionMultiToMulti& f_;
    const Math::RealFunctionMultiDim& errorFunc_;

    double precision_;

    bool updateFile_;
    std::string fileName_;

    int nProcesses_;
    int processId_;

    unsigned long pointsCount_;
    unsigned long newPointsCount_;
    unsigned long newCommunicateCount_;
    unsigned long updateCount_;
    unsigned long minCount_;
    unsigned long communicateCount_;
    unsigned long testSize_;
    unsigned long updateErrorThreshold_;

    int communicateTag_;

    std::vector<void*> updateRequests_;
    std::vector<void*> updateReceiveReq_;
    bool firstUpdateRequested_;

    std::ofstream* logFile_;

    Math::UniformRealGenerator gen_;

    FastApproximator* fa_;
    FastApproximatorError* fast_;

    unsigned long totalCount_, successfulCount_, sameCount_;

    std::vector<std::vector<double> > points_, data_;

    std::vector<double> tempParams_;
    std::vector<double> tempData_;

    std::vector<double> currentParams_;
    std::vector<double> currentData_;

    std::vector<double> communicateBuff_;
    std::vector<std::vector<double> > communicateData_;
    std::vector<std::vector<double> > receiveBuff_;

    std::map<std::vector<double>, unsigned long, PointComp> pointMap_;
};

#endif

