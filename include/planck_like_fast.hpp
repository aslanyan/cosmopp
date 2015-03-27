#ifndef COSMO_PP_PLANCK_LIKE_FAST_HPP
#define COSMO_PP_PLANCK_LIKE_FAST_HPP

#include <fstream>
#include <vector>

#include <planck_like.hpp>
#include <cmb.hpp>
#include <learn_as_you_go.hpp>

/// Fast Planck likelihood class, enhanced by LearnAsYouGo.
class PlanckLikeFast : public Math::LikelihoodFunction
{
public:
    /// Constructor.
    /// \param params A pointer to CosmologicalParameters. This is used to simply set the parameter model and the number of cosmological parameters. The values of the parameters do not matter.
    /// \param useCommander Defines if commander likelihood should be included (true by default).
    /// \param useCamspec Defines if Camspec likelihood should be included (true by default).
    /// \param useLensing Defines if lensing likelihood should be included (true by default).
    /// \param usePolarization Defines if polarization likelihood should be included (false by default).
    /// \param useActSpt Defines if high-l likelihood (ACT and SPT) should be included (false by default).
    /// \param includeTensors Defines if tensor modes should be taken into account during calculations (false by default).
    /// \param kPerDecade The number of points per decade in the k space for the primordial power spectrum calculation.
    /// \param precision The precision of the likelihood. If the estimated error of the approximation is less than this precision then the approximation is used (fast), otherwise the full likelihood will be calculated (slow).
    /// \param minCount The minimum number of points in the training set for the approximation to be acceptable.
    PlanckLikeFast(CosmologicalParams* params, bool useCommander = true, bool useCamspec = true, bool useLensing = true, bool usePolarization = false, bool useActSpt = false, bool includeTensors = false, double kPerDecade = 100, double precision = 0.2, unsigned long minCount = 10000);


    /// Destructor.
    ~PlanckLikeFast();

    /// Calculate the likelihood for a new input point. The calculation will be approximate if the approximation is acceptable.
    /// \param params A vector of the parameters, should always start with the cosmological parameters, followed by camspec extra parameters (if camspec is included), followed by high-l extra parameters (if high l is included).
    /// \param nPar The number of the parameters, used only for checking.
    double calculate(double* params, int nPar) { return doCalculation(params, nPar, false); }

    /// Same as calculate, but enforces exact calculation.
    /// \param params A vector of the parameters, should always start with the cosmological parameters, followed by camspec extra parameters (if camspec is included), followed by high-l extra parameters (if high l is included).
    /// \param nPar The number of the parameters, used only for checking.
    double calculateExact(double* params, int nPar) { return doCalculation(params, nPar, true); }

    /// Set the precision.
    /// \param p The precision of the likelihood. If the estimated error of the approximation is less than this precision then the approximation is used (fast), otherwise the full likelihood will be calculated (slow).
    void setPrecision(double p);

    /// Log the errors for each call of calculate. If this is set then for each call a new row will be added to the log file containing the input parameters followed by the resulting likelihood, the 68.3% and the 95.5% upper bounds of the absolute error, the mean of the error, and the variance of the error.
    /// \param fileNameBase The file name base. If only one process is run then the log file name is simply the base followed by ".txt". If multiple MPI processes are run, each will create a log file with the name "fileNameBase_id.txt", where id is the MPI process ID, i.e. a number between 0 and number of processes - 1.
    void logError(const char* fileNameBase);

private:
    double doCalculation(double* params, int nPar, bool exact);

private:
    CosmologicalParams* cosmoParams_;
    std::vector<double> cosmoParamsVec_;
    int nParams_;

    bool useCommander_;
    bool useCamspec_;
    bool useLensing_;
    bool usePol_;
    bool useActSpt_;

    int lMax_;

    CMB cmb_;
    std::vector<double> lList_;

    PlanckLikelihood like_;
    LearnAsYouGo* layg_;

    std::vector<double> res_;

    void* func_;
    void* errorFunc_;

    std::vector<double> cl_, clTT_;

    std::ofstream outError_;
    bool logError_;
};

#endif

