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

/// Learn as you go approximation class.
/// This class evaluates a given function f, and as it goes it builds a training set. For every new call, it checks whether a quick approximation from the already existing set is acceptable and if so, calculates the approximation. Otherwise the exact value of f is calculated and added to the training set.
/// This class supports parallelism through MPI. Namely, each process will perform the same task, but all of the processes periodically share their training sets with each other, so the training set for each process includes all of the exact calculations by all of the processes.
class LearnAsYouGo
{
public:
    /// Constructor.
    /// \param nIn The dimensionality of the input space.
    /// \param nOut The dimensionality of the output space.
    /// \param f The function to be approximated.
    /// \param errorFunc The function that is used to model the error. errorFunc should take as an input the output of f and return a single real number.
    /// \minCount The minimum size of the training set before approximation can be performed.
    /// \precision The error threshold. This is used to decide whether or not the approximation is acceptable.
    /// \fileName If specified, the training set is continuously saved into this file. Also, if the file exists, the training set will be read in the constructor. So if a file is specified, it will always be read and updated.
    LearnAsYouGo(int nIn, int nData, const Math::RealFunctionMultiToMulti& f, const Math::RealFunctionMultiDim& errorFunc, unsigned long minCount = 10000, double precision = 0.1, const char* fileName = "");

    /// Destructor.
    ~LearnAsYouGo();
    
    /// Evaluate the function for a given point. The approximated result will be returned if the approximation is acceptable, otherwise the exact value will be calculated and returned. If the exact value is calculated it will be added to the training set.
    /// \param x The input point.
    /// \param res The result will be returned here.
    /// \param error1Sigma If specified (i.e. not NULL), the one sigma upper bound of the absolute error probability distribution will be returned here.
    /// \param error2Sigma If specified (i.e. not NULL), the two sigma upper bound of the absolute error probability distribution will be returned here.
    /// \param errorMean If specified (i.e. not NULL), the mean of the error probability distribution will be returned here.
    /// \param errorVar If specified (i.e. not NULL), the variance of the error probability distribution will be returned here.
    void evaluate(const std::vector<double>& x, std::vector<double>* res, double *error1Sigma = NULL, double *error2Sigma = NULL, double *errorMean = NULL, double *errorVar = NULL);

    /// This is similar to evaluate except that it will always calculate the exact value and add to the training set.
    /// \param x The input point.
    /// \param res The result will be returned here.
    void evaluateExact(const std::vector<double>& x, std::vector<double>* res);

    /// Set the precision for the error model.
    /// \p The error threshold. This is used to decide whether or not the approximation is acceptable.
    void setPrecision(double p);

    /// Save into a file.
    /// \param fileName The name of the file.
    void writeIntoFile(const char* fileName) const;
    
    /// Read from a file.
    /// \param fileName The name of the file.
    /// \return If the operation was successful.
    bool readFromFile(const char* fileName);

    /// Set a file to log the progress. Upon every call of evaluate a new row will be added to this log file with the following values: the total number of calls, the number of calls with for which the same input point has been used previously, the number of calls for which the approximation was successful (not counting the cases where the input point was the same as a point in the training set), and the number of calls for which the approximation failed and the exact value of the function was calculated.
    /// \param fileNameBase The file name base. If only one process is run then the log file name is simply the base followed by ".txt". If multiple MPI processes are run, each will create a log file with the name "fileNameBase_id.txt", where id is the MPI process ID, i.e. a number between 0 and number of processes - 1.
    void logIntoFile(const char* fileNameBase);

    /// Get the total number of calls.
    /// \return The total number of calls.
    unsigned long getTotalCount() const { return totalCount_; }

    /// Get the number of successful calls, i.e. calls for which the approximation succeeded.
    /// \return The number of successful calls.
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

