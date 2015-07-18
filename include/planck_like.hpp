#ifndef COSMO_PP_PLANCK_LIKE_HPP
#define COSMO_PP_PLANCK_LIKE_HPP

#include <vector>
#include <string>

#include <likelihood_function.hpp>
#include <cmb.hpp>
#include <random.hpp>

#ifdef COSMO_PLANCK_15

class PlanckLikelihood : public Math::LikelihoodFunction
{
public:
    /// Constructor.
    /// \param lowT Defines if the low-l temperature likelihood should be included.
    /// \param lowP Defines if the low-l polarization likelihood should be included. Note that this can only be inluded if low-l T is included.
    /// \param highT Defines if the high-l temperature likelihood should be included.
    /// \param highP Defines if the high-l polarization likelihood should be included. Note that this can only be inluded if high-l T is included.
    /// \param highLikeLite If true, the "foreground marginalized" version of the high-l likelihoods will be used instead of the full ones.
    /// \param lensingT Defines if the lensing temperature likelihood should be included.
    /// \param lensingP Defines if the lensing polarization likelihood should be included. Note that this can only be inluded if lensing T is included.
    /// \param includeTensors Defines if tensor modes should be taken into account during calculations (false by default).
    /// \param kPerDecade The number of points per decade in the k space for the primordial power spectrum calculation.
    /// \param useOwnCmb If set to true then an instance of the CMB class will be created. In this case the new cosmological parameters can be passed using setCosmoParams. If useOwnCmb is false then the Cl values must be passed through setCls before the likelihoods can be calculated.
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

    /// Set the planck absulute calibration. If not set, the default value of 1 will be used.
    /// \param aPlanck The new absolute calibration value.
    void setAPlanck(double aPlanck);

    /// Set the calibration of the polarization relative to the temperature. If not set, the default value of 1 will be used.
    /// \param aPol The new calibration of the polarization relative to the temperature.
    void setAPol(double aPol);
    
    /// Set the high-l foreground parameters. Should only be called if high-l likelihood is included and the lite version is not used.
    /// \param params A vector containing the parameters. If only temperature is used, there should be 15 parameters in the following order: A_cib_217, cib_index, xi_sz_cib, A_sz, ps_A_100_100, ps_A_143_143, ps_A_143_217, ps_A_217_217, ksz_norm, gal545_A_100, gal545_A_143, gal545_A_143_217, gal545_A_217, calib_100T, calib_217T. If polarization is used as well, there should be 32 parameters in the following order: A_cib_217, cib_index, xi_sz_cib, A_sz, ps_A_100_100, ps_A_143_143, ps_A_143_217, ps_A_217_217, ksz_norm, gal545_A_100, gal545_A_143, gal545_A_143_217, gal545_A_217, galf_EE_A_100, galf_EE_A_100_143, galf_EE_A_100_217, galf_EE_A_143, galf_EE_A_143_217, galf_EE_A_217, galf_EE_index, galf_TE_A_100, galf_TE_A_100_143, galf_TE_A_100_217, galf_TE_A_143, galf_TE_A_143_217, galf_TE_A_217, galf_TE_index, calib_100T, calib_217T, calib_100P, calib_143P, calib_217P. Note that the beam leakage paramters are not included here, and should be set separately using the setBeamLeakageParams function.
    void setHighExtraParams(const std::vector<double>& params);

    /// Set the beam leakage parameters. Should only be called if high-l polarization likelihood is included and the lite version is not used. If not set, the default value of 0 will be used for all of the parameters.
    /// \param params A vector containing the parameters. There should be 60 parameters in the following order: bleak_epsilon_0_0T_0E, bleak_epsilon_1_0T_0E, bleak_epsilon_2_0T_0E, bleak_epsilon_3_0T_0E, bleak_epsilon_4_0T_0E, bleak_epsilon_0_0T_1E, bleak_epsilon_1_0T_1E, bleak_epsilon_2_0T_1E, bleak_epsilon_3_0T_1E, bleak_epsilon_4_0T_1E, bleak_epsilon_0_0T_2E, bleak_epsilon_1_0T_2E, bleak_epsilon_2_0T_2E, bleak_epsilon_3_0T_2E, bleak_epsilon_4_0T_2E, bleak_epsilon_0_1T_1E, bleak_epsilon_1_1T_1E, bleak_epsilon_2_1T_1E, bleak_epsilon_3_1T_1E, bleak_epsilon_4_1T_1E, bleak_epsilon_0_1T_2E, bleak_epsilon_1_1T_2E, bleak_epsilon_2_1T_2E, bleak_epsilon_3_1T_2E, bleak_epsilon_4_1T_2E, bleak_epsilon_0_2T_2E, bleak_epsilon_1_2T_2E, bleak_epsilon_2_2T_2E, bleak_epsilon_3_2T_2E, bleak_epsilon_4_2T_2E, bleak_epsilon_0_0E_0E, bleak_epsilon_1_0E_0E, bleak_epsilon_2_0E_0E, bleak_epsilon_3_0E_0E, bleak_epsilon_4_0E_0E, bleak_epsilon_0_0E_1E, bleak_epsilon_1_0E_1E, bleak_epsilon_2_0E_1E, bleak_epsilon_3_0E_1E, bleak_epsilon_4_0E_1E, bleak_epsilon_0_0E_2E, bleak_epsilon_1_0E_2E, bleak_epsilon_2_0E_2E, bleak_epsilon_3_0E_2E, bleak_epsilon_4_0E_2E, bleak_epsilon_0_1E_1E, bleak_epsilon_1_1E_1E, bleak_epsilon_2_1E_1E, bleak_epsilon_3_1E_1E, bleak_epsilon_4_1E_1E, bleak_epsilon_0_1E_2E, bleak_epsilon_1_1E_2E, bleak_epsilon_2_1E_2E, bleak_epsilon_3_1E_2E, bleak_epsilon_4_1E_2E, bleak_epsilon_0_2E_2E, bleak_epsilon_1_2E_2E, bleak_epsilon_2_2E_2E, bleak_epsilon_3_2E_2E, bleak_epsilon_4_2E_2E.
    void setBeamLeakageParams(const std::vector<double>& params);

    /// Set a prior of the form ksz_norm + 1.6 × A_sz = 9.5±3.0 (the recommended prior by planck). This prior will be included INTO the high-l likelihood. Should only be called if high-l likelihoods are used and not the lite version. If not set, it won't be included by default.
    /// \param szPrior Specifies whether or not the sz prior should be included.
    void setSZPrior(bool szPrior = true);

    /// Calculate the low-l likelihood. Should not be called if low-l likelihoods were not included in the constructor. Must be called after setCosmoParams or setCls.
    /// \return -2ln(likelihood).
    double lowLike();

    /// Calculate the high-l likelihood. Should not be called if high-l likelihoods were not included in the constructor. Must be called after setCosmoParams or setCls.
    /// \return -2ln(likelihood).
    double highLike();

    /// Calculate the lensing likelihood. Should not be called if lensing likelihoods were not included in the constructor. Must be called after setCosmoParams or setCls.
    /// \return -2ln(likelihood).
    double lensingLike();

    /// Calculate all of the likelihoods included in the constructor. Must be called after setCosmoParams or setCls.
    /// \return -2ln(likelihood).
    double likelihood();

    /// Use this parameters to set the model for the calculate function. The number of cosmological parameters will be determined from here, and when calculate is called the cosmological parameters will be assigned to this model.
    /// \param params A pointer to the model parameters. Note that when calculate is called params will be changed to set the new parameters.
    void setModelCosmoParams(CosmologicalParams *params) { modelParams_ = params; modelParams_->getAllParameters(vModel_); }

    /// Calculate the likelihood taking all of the params as an input. This is for the general LikelihoodFunction interface. Can only be called if the model parameters are set by setModelCosmoParams.
    /// \param params A vector of the parameters, should always start with the cosmological parameters, followed by A_planck, followed by the high-l extra parameters (if high-l is included and not the lite version), followed by A_pol (if high-l polarization is included and not the lite version).
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

    bool szPrior_;
    Math::UniformRealGenerator rand_;
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

    Math::UniformRealGenerator rand_;
};

#endif

#endif

