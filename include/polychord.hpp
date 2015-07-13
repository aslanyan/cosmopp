#ifndef COSMO_PP_POLYCHORD_HPP
#define COSMO_PP_POLYCHORD_HPP

#include <vector>
#include <string>

#include <macros.hpp>
#include <likelihood_function.hpp>

/// A PolyChord scanner.

/// This class provides an interface for using PolyChord for scanning parameter regions.
class PolyChord
{
public:
    /// Constructor.
    /// \param nPar The number of parameters.
    /// \param like The likelihood function.
    /// \param nLive The number of live points.
    /// \param fileRoot The root for filenames produced by PolyChord.
    PolyChord(int nPar, Math::LikelihoodFunction& like, int nLive, std::string fileRoot, int nRepeats = 0);

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
    double logLike(double *theta);

private:
    Math::LikelihoodFunction& like_;
    std::vector<double> paramsStarting_, paramsMean_, paramsStd_, paramsBest_, paramsCurrent_;
    std::vector<std::string> paramNames_;
    std::vector<int> priorTypes_;
    std::vector<double> priorMins_, priorMaxs_;
    std::vector<double> paramsFixed_;
    std::vector<bool> isFixed_;
    int nFixed_;
    int n_, nLive_, nRepeats_;
    std::string fileRoot_;

    static bool running_;
};

#endif

