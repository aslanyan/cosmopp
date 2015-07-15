#ifndef COSMO_PP_PLANCK_LIKE_HPP
#define COSMO_PP_PLANCK_LIKE_HPP

#include <vector>
#include <string>

#include <likelihood_function.hpp>
#include <cmb.hpp>

#ifdef COSMO_PLANCK_15

class PlanckLikelihood : public Math::LikelihoodFunction
{
public:
    PlanckLikelihood(bool lowT = true, bool lowP = true, bool highT = true, bool highP = true, bool highLikeLite = true, bool lensingT = true, bool lensingP = true, bool includeTensors = false, double kPerDecade = 100, bool useOwnCmb = true);

    /// Destructor.
    ~PlanckLikelihood();

    /// Set cosmological parameters. Either this function or setCls should be called before likelihood calculation.
    /// \param params The cosmological parameters.
    void setCosmoParams(const CosmologicalParams& params);

    /// Set the Cl values. Either this function or setCosmoParams should be called before likelihood calculation.
    /// \param tt The Cl_TT values.
    /// \param ee The Cl_EE values. Can be NULL if do not want to specify this. Required if polarization likelihood has been initialized.
    /// \param te The Cl_TE values. Can be NULL if do not want to specify this. Required if polarization likelihood has been initialized.
    /// \param pp The Cl_PP values. Can be NULL if do not want to specify this. Required if lensing likelihood has been initialized.
    void setCls(const std::vector<double>* tt, const std::vector<double>* ee = NULL, const std::vector<double>* te = NULL, const std::vector<double> *bb = NULL, const std::vector<double>* pp = NULL);

    void setAPlanck(double aPlanck);
    void setAPol(double aPol);
    void setHighExtraParams(const std::vector<double>& params);
    void setBeamLeakageParams(const std::vector<double>& params);

    double lowLike();
    double highLike();
    double lensingLike();

    /// Calculate all of the likelihoods included in the constructor. Must be called after setCosmoParams or setCls.
    /// \return -2ln(likelihood).
    double likelihood();

    /// Use this parameters to set the model for the calculate function. The number of cosmological parameters will be determined from here, and when calculate is called the cosmological parameters will be assigned to this model.
    /// \param params A pointer to the model parameters. Note that when calculate is called params will be changed to set the new parameters.
    void setModelCosmoParams(CosmologicalParams *params) { modelParams_ = params; modelParams_->getAllParameters(vModel_); }

    /// Calculate the likelihood taking all of the params as an input. This is for the general LikelihoodFunction interface. Can only be called if the model parameters are set by setModelCosmoParams.
    /// \param params A vector of the parameters, should always start with the cosmological parameters, followed by camspec extra parameters (if camspec is included), followed by high-l extra parameters (if high l is included).
    /// \param nPar The number of the parameters, used only for checking.
    /// \return -2ln(likelihood).
    double calculate(double* params, int nPar);

    /// Get the l_max value.
    int getLMax() const { return lMax_; }

private:
    void *low_, *high_, *lens_;
    std::vector<std::string> spectraNames_, lensSpectraNames_;
    int lowLMax_, highLMax_, lensLMax_, lMax_;

    bool lowT_, lowP_, highT_, highP_, highLikeLite_, lensingT_, lensingP_;
    CMB* cmb_;
    double aPlanck_, aPol_;
    std::vector<double> highExtra_, beamExtra_;
    std::vector<double> clTT_, clEE_, clTE_, clBB_, clPP_;
    std::vector<double> prevCosmoParams_;
    std::string prevCosmoParamsName_;
    std::vector<double> currentCosmoParams_;
    double prevLow_, prevLens_;
    bool haveLow_, haveLens_;

    const CosmologicalParams* params_;

    CosmologicalParams* modelParams_;
    std::vector<double> vModel_;

    std::vector<double> input_;
};

#else

/// An interface for calculating the likelihood function for Planck.
class PlanckLikelihood : public Math::LikelihoodFunction
{
public:
    /// Constructor.
    /// \param useCommander Defines if commander likelihood should be included (true by default).
    /// \param useCamspec Defines if Camspec likelihood should be included (true by default).
    /// \param useLensing Defines if lensing likelihood should be included (true by default).
    /// \param usePolarization Defines if polarization likelihood should be included (false by default).
    /// \param useActSpt Defines if high-l likelihood (ACT and SPT) should be included (false by default).
    /// \param includeTensors Defines if tensor modes should be taken into account during calculations (false by default).
    /// \param kPerDecade The number of points per decade in the k space for the primordial power spectrum calculation.
    /// \param useOwnCmb If set to true then an instance of the CMB class will be created. In this case the new cosmological parameters can be passed using setCosmoParams. If useOwnCmb is false then the Cl values must be passed through setCls before the likelihoods can be calculated.
    PlanckLikelihood(bool useCommander = true, bool useCamspec = true, bool useLensing = true, bool usePolarization = false, bool useActSpt = false, bool includeTensors = false, double kPerDecade = 100, bool useOwnCmb = true);

    /// Destructor.
    ~PlanckLikelihood();

    /// Set cosmological parameters. Either this function or setCls should be called before likelihood calculation.
    /// \param params The cosmological parameters.
    void setCosmoParams(const CosmologicalParams& params);

    /// Set the Cl values. Either this function or setCosmoParams should be called before likelihood calculation.
    /// \param tt The Cl_TT values.
    /// \param ee The Cl_EE values. Can be NULL if do not want to specify this. Required if polarization likelihood has been initialized.
    /// \param te The Cl_TE values. Can be NULL if do not want to specify this. Required if polarization likelihood has been initialized.
    /// \param pp The Cl_PP values. Can be NULL if do not want to specify this. Required if lensing likelihood has been initialized.
    void setCls(const std::vector<double>* tt, const std::vector<double>* ee = NULL, const std::vector<double>* te = NULL, const std::vector<double>* pp = NULL);

    /// Set extra parameters needed for Camspec. Needs to be called before likelihood calculation. Should not be called if Camspec likelihood is not included.
    void setCamspecExtraParams(double A_ps_100, double A_ps_143, double A_ps_217, double A_cib_143, double A_cib_217, double A_sz, double r_ps, double r_cib, double n_Dl_cib, double cal_100, double cal_217, double xi_sz_cib, double A_ksz, double Bm_1_1);

    /// Set extra parameters needed for high-l (ACT and SPT). Needs to be called before likelihood calculation. Should not be called if high-l likelihood is not included.
    void setActSptExtraParams(double A_sz, double A_ksz, double xi_sz_cib, double a_ps_act_148, double a_ps_act_217, double a_ps_spt_95, double a_ps_spt_150, double a_ps_spt_220, double A_cib_143, double A_cib_217, double n_Dl_cib, double r_ps_spt_95x150, double r_ps_spt_95x220, double r_ps_150x220, double r_cib, double a_gs, double a_ge, double cal_acts_148, double cal_acts_217, double cal_acte_148, double cal_acte_217, double cal_spt_95, double cal_spt_150, double cal_spt_220);

    /// Calculate commander likelihood. Should not be called if commander likelihood was not included in the constructor. Must be called after setCosmoParams or setCls.
    /// \return -2ln(likelihood).
    double commanderLike();

    /// Calculate Camspec likelihood. Should not be called if Camspec likelihood was not included in the constructor. Must be called after setCosmoParams or setCls.
    /// \return -2ln(likelihood).
    double camspecLike();

    /// Calculate polarization likelihood. Should not be called if polarization likelihood was not included in the constructor. Must be called after setCosmoParams or setCls.
    /// \return -2ln(likelihood).
    double polLike();

    /// Calculate lensing likelihood. Should not be called if lensing likelihood was not included in the constructor. Must be called after setCosmoParams or setCls.
    /// \return -2ln(likelihood).
    double lensingLike();

    /// Calculate high-l (ACT and SPT) likelihood. Should not be called if high-l likelihood was not included in the constructor. Must be called after setCosmoParams or setCls.
    /// \return -2ln(likelihood).
    double actSptLike();

    /// Calculate all of the likelihoods included in the constructor. Must be called after setCosmoParams or setCls.
    /// \return -2ln(likelihood).
    double likelihood();

    /// Use this parameters to set the model for the calculate function. The number of cosmological parameters will be determined from here, and when calculate is called the cosmological parameters will be assigned to this model.
    /// \param params A pointer to the model parameters. Note that when calculate is called params will be changed to set the new parameters.
    void setModelCosmoParams(CosmologicalParams *params) { modelParams_ = params; modelParams_->getAllParameters(vModel_); }

    /// Calculate the likelihood taking all of the params as an input. This is for the general LikelihoodFunction interface. Can only be called if the model parameters are set by setModelCosmoParams.
    /// \param params A vector of the parameters, should always start with the cosmological parameters, followed by camspec extra parameters (if camspec is included), followed by high-l extra parameters (if high l is included).
    /// \param nPar The number of the parameters, used only for checking.
    /// \return -2ln(likelihood).
    double calculate(double* params, int nPar);

    /// Get the l_max value.
    int getLMax() const { return lMax_; }

private:
    void* commander_, *camspec_, *pol_, *actspt_;
    void* lens_;
    int commanderLMax_, camspecLMax_, lensingLMax_, polLMax_, actSptLMax_, lMax_;
    CMB* cmb_;
    std::vector<std::string> spectraNames_;
    std::vector<double> camspecExtra_, actSptExtra_;
    std::vector<double> clTT_, clEE_, clTE_, clPP_;
    std::vector<double> prevCosmoParams_;
    std::string prevCosmoParamsName_;
    std::vector<double> currentCosmoParams_;
    double prevCommander_, prevPol_, prevLens_;
    bool haveCommander_, havePol_, haveLens_;

    const CosmologicalParams* params_;

    CosmologicalParams* modelParams_;
    std::vector<double> vModel_;
};

#endif

#endif

