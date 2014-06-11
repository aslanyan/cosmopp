#ifndef COSMO_PP_WMAP9_LIKE_HPP
#define COSMO_PP_WMAP9_LIKE_HPP

#include <vector>
#include <string>

#include <macros.hpp>
#include <likelihood_function.hpp>
#include <cmb.hpp>

/// An interface for calculating the likelihood function for WMAP9.
class WMAP9Likelihood : public Math::LikelihoodFunction
{
public:
    /// Constructor.

    /// Initializes the WMAP9 likelihood code. Note that this code can be initialized only once during the lifetime of the program, so there should be only one WMAP9Likelihood object during the whole run of the program.
    /// \param useLowlT Specifies if low l temperature likelihood should be included.
    /// \param useHighlT Specifies if high l temperature likelihood should be included.
    /// \param useLowlP Specifies if low l polarization likelihood should be included.
    /// \param useHighlP Specifies if high l polarization likelihood should be included.
    /// \param useGibbs Specifies if the Gibbs likelihood should be used for low l temperature calculation (relevant only if it's included).
    /// \param useTTBeam Specifies if the beam and point source likelihood correction should be calculated for temperature.
    WMAP9Likelihood(bool useLowlT = true, bool useHighlT = true, bool useLowlP = true, bool useHighlP = true, bool useGibbs = true, bool useTTBeam = true);

    /// Destructor.
    ~WMAP9Likelihood();

    /// Set cosmological parameters. Needs to be called before calculateCls.
    /// \param params The cosmological parameters.
    void setCosmoParams(const CosmologicalParams& params);

    /// Calculate the Cl values. Should be called after setCosmoParams but before any likelihood calculation.
    void calculateCls();

    /// Calculate all of the likelihoods included in the constructor. Must be called after calculateCls.
    /// \return -2ln(likelihood).
    double likelihood();

    /// Calculate the likelihood taking all of the params as an input. This is for the general LikelihoodFunction interface.
    /// \param params A vector of the parameters, this implementation assumes the 6 LCDM parameters in the order ombh2, omch2, h, tau, ns, As (assumed pivot is 0.05Mpc^-1).
    /// \param nPar The number of the parameters, used only for checking. In this implementation this should always be 6.
    /// \return -2ln(likelihood).
    double calculate(double* params, int nPar);

    /// Get the low-l TT likelihood. Should be called after the likelihood function is called to calculate all of the likelihoods. Can be called only if useLowlT is true in the constructor.
    /// \return -2ln(likelihood).
    inline double lowlTLike() const { check(useLowlT_, "low l T likelihood not initialized in constructor"); return lowlTChi2_ + lowlTDet_; }

    /// Get the low-l TT chi squared. Should be called after the likelihood function is called to calculate all of the likelihoods. Can be called only if useLowlT is true in the constructor.
    /// \return chi^2.
    inline double lowlTChi2() const { check(useLowlT_, "low l T likelihood not initialized in constructor"); return lowlTChi2_; }

    /// Get the low-l TT covariance matrix determinant. Should be called after the likelihood function is called to calculate all of the likelihoods. Can be called only if useLowlT is true in the constructor.
    /// \return ln(det).
    inline double lowlTDet() const { check(useLowlT_, "low l T likelihood not initialized in constructor"); return lowlTDet_; }

    /// Get the high-l TT likelihood. Should be called after the likelihood function is called to calculate all of the likelihoods. Can be called only if useHighlT is true in the constructor.
    /// \return -2ln(likelihood).
    inline double highlTLike() const { check(useHighlT_, "high l T likelihood not initialized in constructor"); return highlT_; }

    /// Get the low-l Polarization likelihood. Should be called after the likelihood function is called to calculate all of the likelihoods. Can be called only if useLowlP is true in the constructor.
    /// \return -2ln(likelihood).
    inline double lowlPLike() const { check(useLowlP_, "low l P likelihood not initialized in constructor"); return lowlPChi2_ + lowlPDet_; }

    /// Get the low-l Polarization chi squared. Should be called after the likelihood function is called to calculate all of the likelihoods. Can be called only if useLowlP is true in the constructor.
    /// \return chi^2.
    inline double lowlPChi2() const { check(useLowlP_, "low l P likelihood not initialized in constructor"); return lowlPChi2_; }

    /// Get the low-l Polarization covariance matrix determinant. Should be called after the likelihood function is called to calculate all of the likelihoods. Can be called only if useLowlP is true in the constructor.
    /// \return ln(det).
    inline double lowlPDet() const { check(useLowlP_, "low l P likelihood not initialized in constructor"); return lowlPDet_; }

    /// Get the high-l TE likelihood. Should be called after the likelihood function is called to calculate all of the likelihoods. Can be called only if useHighlP is true in the constructor.
    /// \return -2ln(likelihood).
    inline double TELike() const { check(useHighlP_, "high l P likelihood not initialized in constructor"); return teChi2_ + teDet_; }

    /// Get the high-l TE chi squared. Should be called after the likelihood function is called to calculate all of the likelihoods. Can be called only if useHighlP is true in the constructor.
    /// \return chi^2.
    inline double TEChi2() const { check(useHighlP_, "high l P likelihood not initialized in constructor"); return teChi2_; }

    /// Get the high-l TE covariance matrix determinant. Should be called after the likelihood function is called to calculate all of the likelihoods. Can be called only if useHighlP is true in the constructor.
    /// \return ln(det).
    inline double TEDet() const { check(useHighlP_, "high l P likelihood not initialized in constructor"); return teDet_; }

    /// Get the high-l TB likelihood. Should be called after the likelihood function is called to calculate all of the likelihoods. Can be called only if useHighlP is true in the constructor.
    /// \return -2ln(likelihood).
    inline double TBLike() const { check(useHighlP_, "high l P likelihood not initialized in constructor"); return tbChi2_ + tbDet_; }

    /// Get the high-l TB chi squared. Should be called after the likelihood function is called to calculate all of the likelihoods. Can be called only if useHighlP is true in the constructor.
    /// \return chi^2.
    inline double TBChi2() const { check(useHighlP_, "high l P likelihood not initialized in constructor"); return tbChi2_; }

    /// Get the high-l TB covariance matrix determinant. Should be called after the likelihood function is called to calculate all of the likelihoods. Can be called only if useHighlP is true in the constructor.
    /// \return ln(det).
    inline double TBDet() const { check(useHighlP_, "high l P likelihood not initialized in constructor"); return tbDet_; }

    /// Get the TT beam and point source correction likelihood. Should be called after the likelihood function is called to calculate all of the likelihoods. Can be called only if useTTBeam is true in the constructor.
    /// \return -2ln(likelihood).
    inline double TTBeamLike() const { check(useTTBeam_, "TT beam likelihood not initialized in constructor"); return ttBeam_; }

private:
    CMB cmb_;
    std::vector<double> clTT_, clEE_, clTE_, clTB_, clEB_, clBB_;
    std::vector<double> like_;
    bool useLowlT_, useHighlT_, useLowlP_, useHighlP_, useGibbs_, useTTBeam_;
    double lowlTChi2_, lowlTDet_, highlT_, teChi2_, teDet_, lowlPChi2_, lowlPDet_, tbChi2_, tbDet_, ttBeam_;
    static bool initialized_;
};

#endif

