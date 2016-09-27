#ifndef COSMO_PP_PLANCKSZ_LIKE_HPP
#define COSMO_PP_PLANCKSZ_LIKE_HPP

#include <string>
#include <fstream>

#include <macros.hpp>
#include <cmb.hpp>
#include <likelihood_function.hpp>

class PlanckSZLikelihood : public Math::LikelihoodFunction
{
public:
    PlanckSZLikelihood(CMB& cmb, bool initializeCosmoAtEachStep = true, bool fixedBias = false) : initializeCosmoAtEachStep_(initializeCosmoAtEachStep), fixedBias_(fixedBias), cmb_(&cmb)
    {
    }

    ~PlanckSZLikelihood() {}

    void setCosmoParams(const CosmologicalParams& params)
    {
        params_ = &params;
        if(initializeCosmoAtEachStep_)
            cmb_->initialize(params, true, false, false, true);
    }

    double likelihood()
    {
        double mean = 0.764;
        double stddev = 0.025;
        if(fixedBias_)
        {
            mean = 0.78;
            stddev = 0.01;
        }
        double sigma8 = cmb_->sigma8();
        double OmM = params_->getOmM();
        double LnLike = pow(((sigma8*pow(OmM/0.27, 0.3) - mean)/stddev), 2);
        return LnLike;
    }

    void setModelCosmoParams(CosmologicalParams *params)
    {
        // Sets modelParams_ to initial cosmological parameters.
        modelParams_ = params;
        // Sets vModel_ to a vector containing the parameters in modelParams_.
        modelParams_->getAllParameters(vModel_);
    }

    double calculate(double* params, int nPar)
    {
        // Check to see that modelParams_ is set. This is set by calling setModelCosmoParams().
        check(modelParams_, "model params must be set before calling this function");
        // Check to see that vModel_ is not empty
        check(!vModel_.empty(), "");
        const int nModel = vModel_.size();
        
        // Check that number of parameters passed in is same as number of parameters in nModel
        check(nPar == nModel, "");
    
        bool modelHasChanged = false;
        // Set all the parameters in vModel_ to the values in params
        for(int i = 0; i < nModel; ++i)
        {
            if(vModel_[i] != params[i])
            {
                modelHasChanged = true;
                vModel_[i] = params[i];
            }
        }

        // Sets the parameters in modelParams_ to the values in vModel_
        if(modelHasChanged)
        {
            double badLike = 0;
            const bool success = modelParams_->setAllParameters(vModel_, &badLike);
            if(!success)
            {
                check(badLike >= 0, "");
                return 1e10 + badLike;
            }

            check(badLike == 0, "");

            // Set the cosmological parameters to modelParams_.
            setCosmoParams(*modelParams_);
        }
    
        return likelihood();
    }

private:
    CMB* cmb_;
    const bool initializeCosmoAtEachStep_;

    const bool fixedBias_;

    const CosmologicalParams* params_;
    
    CosmologicalParams* modelParams_;
    std::vector<double> vModel_; // modelParams_ as a vector
};

#endif

