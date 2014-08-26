#ifndef COSMO_PP_PLANCK_LIKE_HPP
#define COSMO_PP_PLANCK_LIKE_HPP

#include <vector>
#include <string>

#include <likelihood_function.hpp>
#include <cmb.hpp>

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
    PlanckLikelihood(bool useCommander = true, bool useCamspec = true, bool useLensing = true, bool usePolarization = false, bool useActSpt = false, bool includeTensors = false, double kPerDecade = 100);

    /// Destructor.
    ~PlanckLikelihood();

    /// Sets a reference to a CMB object which will be used to get the power spectra. If this is set then setCosmoParams should not be used any more (will not affect the likelihood).
    /// \param cmb A reference to the CMB object.
    void setCMB(CMB& cmb) { useCMB_ = &cmb; }

    /// Set cosmological parameters. Needs to be called before calculateCls.
    /// \param params The cosmological parameters.
    void setCosmoParams(const CosmologicalParams& params);

    /// Set extra parameters needed for Camspec. Needs to be called before likelihood calculation. Should not be called if Camspec likelihood is not included.
    void setCamspecExtraParams(double A_ps_100, double A_ps_143, double A_ps_217, double A_cib_143, double A_cib_217, double A_sz, double r_ps, double r_cib, double n_Dl_cib, double cal_100, double cal_217, double xi_sz_cib, double A_ksz, double Bm_1_1);

    /// Set extra parameters needed for high-l (ACT and SPT). Needs to be called before likelihood calculation. Should not be called if high-l likelihood is not included.
    void setActSptExtraParams(double A_sz, double A_ksz, double xi_sz_cib, double a_ps_act_148, double a_ps_act_217, double a_ps_spt_95, double a_ps_spt_150, double a_ps_spt_220, double A_cib_143, double A_cib_217, double n_Dl_cib, double r_ps_spt_95x150, double r_ps_spt_95x220, double r_ps_150x220, double r_cib, double a_gs, double a_ge, double cal_acts_148, double cal_acts_217, double cal_acte_148, double cal_acte_217, double cal_spt_95, double cal_spt_150, double cal_spt_220);

    /// Calculate the Cl values. Should be called after setCosmoParams but before any likelihood calculation.
    void calculateCls();

    /// Calculate commander likelihood. Should not be called if commander likelihood was not included in the constructor. Must be called after calculateCls.
    /// \return -2ln(likelihood).
    double commanderLike();

    /// Calculate Camspec likelihood. Should not be called if Camspec likelihood was not included in the constructor. Must be called after calculateCls.
    /// \return -2ln(likelihood).
    double camspecLike();

    /// Calculate polarization likelihood. Should not be called if polarization likelihood was not included in the constructor. Must be called after calculateCls.
    /// \return -2ln(likelihood).
    double polLike();

    /// Calculate lensing likelihood. Should not be called if lensing likelihood was not included in the constructor. Must be called after calculateCls.
    /// \return -2ln(likelihood).
    double lensingLike();

    /// Calculate high-l (ACT and SPT) likelihood. Should not be called if high-l likelihood was not included in the constructor. Must be called after calculateCls.
    /// \return -2ln(likelihood).
    double actSptLike();

    /// Calculate all of the likelihoods included in the constructor. Must be called after calculateCls.
    /// \return -2ln(likelihood).
    double likelihood();

    /// Calculate the likelihood taking all of the params as an input. This is for the general LikelihoodFunction interface.
    /// \param params A vector of the parameters, should always start with 6 LCDM parameters ombh2, omch2, h, tau, ns, As (assumed pivot is 0.05Mpc^-1), followed by camspec extra parameters (if camspec is included), followed by high-l extra parameters (if high l is included).
    /// \param nPar The number of the parameters, used only for checking.
    /// \return -2ln(likelihood).
    double calculate(double* params, int nPar);

private:
    void* commander_, *camspec_, *pol_, *actspt_;
    void* lens_;
    int commanderLMax_, camspecLMax_, lensingLMax_, polLMax_, actSptLMax_, lMax_;
    CMB cmb_;
    CMB* useCMB_;
    std::vector<std::string> spectraNames_;
    std::vector<double> camspecExtra_, actSptExtra_;
    std::vector<double> clTT_, clEE_, clTE_, clPP_;
    std::vector<double> cosmoParams_;
    bool prevCosmoCalculated_;
};

#endif

