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
    /// \param speed The speed grade of the parameter. Should be 1 or higher. Default value is 1.
    void setParam(int i, const std::string& name, double min, double max, int speed = 1);

    /// Define a given parameter to have a log uniform prior. One of the parameter setting functions must be called for each parameter before the run.
    /// \param i The index of the parameter, 0 <= i < number of parameters.
    /// \param name The name of the parameter.
    /// \param min The minimum value of the parameter (the lower bound for the prior).
    /// \param max The maximum value of the parameter (the upper bound for the prior).
    /// \param speed The speed grade of the parameter. Should be 1 or higher. Default value is 1.
    void setParamLogUniform(int i, const std::string& name, double min, double max, int speed = 1);

    /// Define a given parameter to have a gaussian prior. One of the parameter setting functions must be called for each parameter before the run.
    /// \param i The index of the parameter, 0 <= i < number of parameters.
    /// \param name The name of the parameter.
    /// \param mean The mean of the prior
    /// \param sigma The sigma of the prior
    /// \param speed The speed grade of the parameter. Should be 1 or higher. Default value is 1.
    void setParamGauss(int i, const std::string& name, double mean, double sigma, int speed = 1);

    /// Define a given parameter to have a sorted uniform prior.
    /// This means that this parameter is a part of a group of parameters which have a joint prior.
    /// They are distributed uniformly but with the constraint that they are sorted. A parameter with a smaller index will always be smaller than a parameter with a larger index.
    /// One of the parameter setting functions must be called for each parameter before the run.
    /// \param i The index of the parameter, 0 <= i < number of parameters.
    /// \param name The name of the parameter.
    /// \param min The minimum value of the parameter (the lower bound for the prior).
    /// \param max The maximum value of the parameter (the upper bound for the prior).
    /// \param block This is an index that determines the block of parameters with the sorted uniform prior. All the parameters to be a part of the same group of sorted parameters must have the same block. The block index must be non-negative.
    /// \param speed The speed grade of the parameter. Should be 1 or higher. Default value is 1.
    void setParamSortedUniform(int i, const std::string& name, double min, double max, int block, int speed = 1);

    /// Define a given parameter to be fixed to a given value.
    /// \param i The index of the parameter, 0 <= i < number of parameters.
    /// \param name The name of the parameter.
    /// \param val The fixed value of the parameter.
    void setParamFixed(int i, const std::string& name, double val);

    /// Get the name of a parameter.
    /// \param i The index of the parameter.
    /// \return The name of the parameter.
    const std::string& getParamName(int i) const { check(i >= 0 && i < n_, "invalid index " << i); return paramNames_[i]; }

    /// Set the parameter hierarchy.
    /// \parm fracs A vector that contains the fractions of time spent in each hierarchy level (the sum should be one).
    void setParameterHierarchy(const std::vector<double>& fracs);

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
    std::vector<int> priorBlocks_;
    std::vector<double> priorMins_, priorMaxs_;
    std::vector<int> speeds_;
    std::vector<double> paramsFixed_;
    std::vector<bool> isFixed_;
    int nFixed_;
    int n_, nLive_, nRepeats_;
    std::string fileRoot_;

    static bool running_;

    std::vector<double> fracs_;
};

#endif

