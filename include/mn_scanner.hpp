#ifndef COSMO_PP_MN_SCANNER_HPP
#define COSMO_PP_MN_SCANNER_HPP

#include <vector>
#include <string>

#include <likelihood_function.hpp>
#include <table_function.hpp>

/// A Multinest scanner.

/// This class provides an interface for using Multinest for scanning parameter regions.
class MnScanner
{
public:
    /// Constructor.
    /// \param nPar The number of parameters.
    /// \param like The likelihood function.
    /// \param nLive The number of live points.
    /// \param fileRoot The root for filenames produced by Multinest.
    /// \param accurateEvidence Defines if accurate Bayesian evidence calculation is required (true by default). If this is not required the run will be faster.
    MnScanner(int nPar, Math::LikelihoodFunction& like, int nLive, std::string fileRoot, bool accurateEvidence = true);

    /// Define a given parameter to have a uniform prior. One of the parameter setting functions must be called for each parameter before the run.
    /// \param i The index of the parameter, 0 <= i < number of parameters.
    /// \param name The name of the parameter.
    /// \param min The minimum value of the parameter (the lower bound for the prior).
    /// \param max The maximum value of the parameter (the upper bound for the prior).
    void setParam(int i, const std::string& name, double min, double max);

    /// Define a given parameter to have a gaussian prior. One of the parameter setting functions must be called for each parameter before the run.
    /// \param i The index of the parameter, 0 <= i < number of parameters.
    /// \param name The name of the parameter.
    /// \param mean The mean of the prior
    /// \param sigma The sigma of the prior
    void setParamGauss(int i, const std::string& name, double mean, double sigma);

    /// Define a given parameter to have a general prior. One of the parameter setting functions must be called for each parameter before the run.
    /// \param i The index of the parameter, 0 <= i < number of parameters.
    /// \param name The name of the parameter.
    /// \param min The minimum value of the parameter (the lower bound for the prior).
    /// \param max The maximum value of the parameter (the upper bound for the prior).
    /// \param distrib The prior function. This function should be valid between min and max.
    void setParamGeneral(int i, const std::string& name, double min, double max, const Math::RealFunction& distrib);

    /// Define a given parameter to be fixed to a given value.
    /// \param i The index of the parameter, 0 <= i < number of parameters.
    /// \param name The name of the parameter.
    /// \param val The fixed value of the parameter.
    void setParamFixed(int i, const std::string& name, double val);

    /// Get the name of a parameter.
    /// \param i The index of the parameter.
    /// \return The name of the parameter.
    const std::string& getParamName(int i) const { check(i >= 0 && i < n_, "invalid index " << i); return paramNames_[i]; }

    /// Run the scan. Should be called after all of the other necessary functions have been called to set all of the necessary settings. The resulting chain is written in the file (fileRoot).txt. The first column is the number of repetitions of the element, the second column is -2ln(likelihood), the following columns are the values of all of the parameters.
    /// \param resume Resume from previous job or not (true by default).
    void run(bool resume = true);

public:
    void logLike(double *Cube, int &ndim, int &npars, double &lnew);
    void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &logZerr);

private:
    void dumpInfo(const char* error);

private:
    Math::LikelihoodFunction& like_;
    std::vector<double> paramsStarting_, paramsMean_, paramsStd_, paramsBest_, paramsCurrent_;
    std::vector<std::string> paramNames_;
    std::vector<Math::TableFunction<double, double> > paramPriors_;
    std::vector<double> paramsFixed_;
    int nFixed_;
    int n_, nLive_;
    std::string fileRoot_;
    bool accurateEvidence_;
};

#endif

