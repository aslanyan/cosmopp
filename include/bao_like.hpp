#ifndef COSMO_PP_BAO_LIKE_HPP
#define COSMO_PP_BAO_LIKE_HPP

#include <string>
#include <fstream>

#include <macros.hpp>
#include <likelihood_function.hpp>
#include <cmb.hpp>

/// BAO Likelihood (written by Nicolas Canac)
class BAOLikelihood : public Math::CosmoLikelihood
{
public:
    /// Constructor.
    /// \param cmb A reference to a CMB object to be used.
    /// \param initializeCMBAtEachStep Determines if the cmb object should be initialized every time setCosmoParams is called. Note that the BAO likelihood may be used in combination with other likelihoods, such as Planck, for some kind of parameter space scanning (MCMC, MultiNest). In this case the slowest part is the initialization of CMB. It will make sense to have one CMB object serve multiple likelihoods, so the initialize of CMB will be called outside. If this is the case then this parameter needs to be set to false. The assumption then will be that the CMB object will be initialized with the same cosmological parameters as setCosmoParams.
    /// \param sixDF Determines if the 6DF likelihood should be included.
    /// \param lowZ Determines if the BOSS lowZ DR 10-11 likelihood should be included.
    /// \param cmass Determines if the BOSS CMASS DR 10-11 likelihood should be included.
    /// \param mgs Determines if the MGS DR7 likelihood should be included.
    BAOLikelihood(CMB& cmb, bool initializeCMBAtEachStep = true, bool sixDF = true, bool lowZ = true, bool cmass = true, bool mgs = true) : initializeCMBAtEachStep_(initializeCMBAtEachStep), cmb_(&cmb), sixDF_(sixDF), lowZ_(lowZ), cmass_(cmass), mgs_(mgs)
    {
    }

    /// Destructor.
    ~BAOLikelihood() {}

    /// Set the cosmological parameters. This sets the parameters for which the bao likelihoods will be subsequently calculated.
    /// The initializeCMBAtEachStep parameter of the constructor is important for this function. If that parameter was set to true
    /// then the cmb object will be initialized with these parameters. However, if it was set to false then it is expected that the cmb
    /// object will be initialized separately outside this function with the same parameters.
    /// \param params The cosmological parameters object.
    void setCosmoParams(const CosmologicalParams& params)
    {
        params_ = &params;
        if(initializeCMBAtEachStep_)
            cmb_->initialize(params, true, false, false, false);
    }

    /// 6DF likelihood
    /// \return -2ln(likelihood) for 6DF
    double SixDF_Likelihood()
    {
        double da, dr, dv, rsdrag, theo;
        double z, value, error;

        // 6DF
        z = 0.106;
        value = 0.327; // rs/D_V
        error = 0.015;

        da = cmb_->getAngularDistance(z);
        dr = z / cmb_->getHubble(z);
        dv = pow(da*da*(1 + z)*(1 + z)*dr,1.0/3.0);
        rsdrag = cmb_->getrsdrag();
        theo = rsdrag/dv;
        return pow((theo - value)/error, 2.0);
    }

    /// BOSS LowZ DR10&11 likelihood (Anderson et al. 1312.4877)
    /// \return -2ln(likelihood) for BOSS LowZ
    double LowZ_DR10_11_Likelihood()
    {
        double da, dr, dv, rsdrag, theo;
        double z, value, error;

        z = 0.32;
        value = 8.47; // D_V/rs
        error = 0.17;

        da = cmb_->getAngularDistance(z);
        dr = z / cmb_->getHubble(z);
        dv = pow(da*da*(1 + z)*(1 + z)*dr,1.0/3.0);
        rsdrag = cmb_->getrsdrag();
        theo = dv/rsdrag;
        return pow((theo - value)/error, 2.0);
    }

    /// BOSS CMASS DR10&11 likelihood (Anderson et al. 1312.4877)
    /// \return -2ln(likelihood) for BOSS CMASS
    double CMASS_DR10_11_Likelihood()
    {
        double da, dr, dv, rsdrag, theo;
        double z, value, error;

        z = 0.57;
        value = 13.77; // D_V/rs
        error = 0.13;

        da = cmb_->getAngularDistance(z);
        dr = z / cmb_->getHubble(z);
        dv = pow(da*da*(1 + z)*(1 + z)*dr,1.0/3.0);
        rsdrag = cmb_->getrsdrag();
        theo = dv/rsdrag;
        return pow((theo - value)/error, 2.0);
    }

    /// MGS SDSS DR7 likelihood (Ross et al. 1409.3242v1)
    /// \return -2ln(likelihood) for MGS
    double MGS_DR7_Likelihood()
    {
        double da, dr, dv, rsdrag, theo;
        double z, value, error;

        z = 0.15;
        value = 4.47; // D_V/rs
        error = 0.16;

        da = cmb_->getAngularDistance(z);
        dr = z / cmb_->getHubble(z);
        dv = pow(da*da*(1 + z)*(1 + z)*dr,1.0/3.0);
        rsdrag = cmb_->getrsdrag();
        theo = dv/rsdrag;
        return pow((theo - value)/error, 2.0);
    }

    /// Full likelihood
    /// \return -2ln(likelihood) where likelihood is the product of a few different likelihoods. The parameters of the constructor determine which likelihoods will be included.
    double likelihood()
    {
        double chi2 = 0;
        if(sixDF_) chi2 += SixDF_Likelihood();
        if(lowZ_) chi2 += LowZ_DR10_11_Likelihood();
        if(cmass_) chi2 += CMASS_DR10_11_Likelihood();
        if(mgs_) chi2 += MGS_DR7_Likelihood();
        return chi2;
    }

    double SDSSLikelihood()
    {
        double c = 2.99792458e8;
        // From data set in http://arxiv.org/abs/0907.1660
        double z1 = 0.2, z2 = 0.35;
        double rstodvz1 = 0.190533, rstodvz2 = 0.109715;
        double invcov[2][2];
        invcov[0][0] = 30124.1;
        invcov[0][1] = -17226.0;
        invcov[1][0] = invcov[0][1];
        invcov[1][1] = 86976.6;

        double OmK, OmLambda, OmM, h, H0, w, rsdrag;
        double dv1theory, dv2theory, hz1, hz2, DAz1, DAz2;
        double rstodvz1theorydelta, rstodvz2theorydelta, LnLike;

        h = params_->getH();
        H0 = 100.0*h;
        OmK = params_->getOmK();
        OmLambda = params_->getOmLambda();
        OmM = params_->getOmM();
        w = -1.0;
        rsdrag = cmb_->getrsdrag();
        DAz1 = cmb_->getAngularDistance(z1);
        DAz2 = cmb_->getAngularDistance(z2);

        hz1 = sqrt( OmM*pow(1.0+z1,3.0) + OmK*pow(1.0+z1,2.0) + OmLambda*pow(1.0+z1,3.0*(1.0+w)) );
        hz2 = sqrt( OmM*pow(1.0+z2,3.0) + OmK*pow(1.0+z2,2.0) + OmLambda*pow(1.0+z2,3.0*(1.0+w)) );
        dv1theory = pow(((1.0+z1)*DAz1),2.0)*c*z1/H0/hz1/1000.0;
        dv2theory = pow(((1.0+z2)*DAz2),2.0)*c*z2/H0/hz2/1000.0;
        dv1theory = pow(dv1theory,(1.0/3.0));
        dv2theory = pow(dv2theory,(1.0/3.0));

        rstodvz1theorydelta = rsdrag/dv1theory - rstodvz1;
        rstodvz2theorydelta = rsdrag/dv2theory - rstodvz2;

        LnLike = 0.5*((rstodvz1theorydelta) * invcov[0][0] * (rstodvz1theorydelta)
                + 2.0 * (rstodvz1theorydelta) * invcov[0][1] * (rstodvz2theorydelta)
                + (rstodvz2theorydelta) * invcov[1][1] * (rstodvz2theorydelta));
    
        return LnLike;
    }

    /// Set the model cosmological parameters. This function must be called before the calculate function.
    /// The model cosmological parameters will be assumed to be the same model as the parameters passed into calculate.
    void setModelCosmoParams(CosmologicalParams *params)
    {
        // Sets modelParams_ to initial cosmological parameters.
        modelParams_ = params;
        // Sets vModel_ to a vector containing the parameters in modelParams_.
        modelParams_->getAllParameters(vModel_);
    }

    /// Calculate the likelihood. Note that setModelCosmoParams must be called before this function.
    /// \param params A pointer to the array of parameters.
    /// \param nPar The number of parameters.
    /// \return -2ln(likelihood).
    double calculate(double* params, int nPar)
    {
        // Check to see that modelParams_ is set. This is set by calling setModelCosmoParams().
        check(modelParams_, "model params must be set before calling this function");
        // Check to see that vModel_ is not empty
        check(!vModel_.empty(), "");
        const int nModel = vModel_.size();
        
        // Check that number of parameters passed in is same as number of parameters in nModel
        check(nPar == nModel, "");
    
        // Set all the parameters in vModel_ to the values in params
        for(int i = 0; i < nModel; ++i)
            vModel_[i] = params[i];
    
        // Sets the parameters in modelParams_ to the values in vModel_
        // NEED TO DO: Write setAllParameters for my cosmological parameters class
        modelParams_->setAllParameters(vModel_);
        // Set the cosmological parameters to modelParams_.
        setCosmoParams(*modelParams_);
    
        return likelihood();
    }

private:
    CMB* cmb_;
    const bool initializeCMBAtEachStep_;
    const bool sixDF_, lowZ_, cmass_, mgs_;

    const CosmologicalParams* params_; // Standard cosmological parameters
    
    CosmologicalParams* modelParams_; // Cosmological parameters in test models
    std::vector<double> vModel_; // modelParams_ as a vector
};

#endif
