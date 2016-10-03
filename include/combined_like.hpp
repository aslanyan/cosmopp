#ifndef COSMO_PP_COMBINED_LIKE_HPP
#define COSMO_PP_COMBINED_LIKE_HPP

#include <string>
#include <fstream>
#include <memory>

#include <macros.hpp>
#include <planck_like.hpp>
#include <wmap9_like.hpp>
#include <bao_like.hpp>
#include <numerics.hpp>
#include <wigglez_like.hpp>
#include <plancksz_like.hpp>
#include <h0_like.hpp>

/// Combined Likelihood (written by Nicolas Canac)
class CombinedLikelihood : public Math::CosmoLikelihood
{
public:
    CombinedLikelihood(std::string datapath, bool usePlanck, bool useBAO, bool useWiggleZ, bool useSZ = false, bool useH0 = false, bool useHighP = false, bool planckLikeLite = true) : usePlanck_(usePlanck), useBAO_(useBAO), useWiggleZ_(useWiggleZ), useSZ_(useSZ), useH0_(useH0), useHighP_(useHighP), planckLikeLite_(planckLikeLite), planckExtraParams_(planckLikeLite ? 0 : 15 + (useHighP ? 17 : 0))
    {
        cmb_.preInitialize(3500, false, true, false, 0, 100, 1e-6, 1.0);
        //check(!(usePlanck_ && useWMAP_), "Both Planck and WMAP likelihoods should not be used at the same time.");
        //check(!(useBAO_ && useLRG_), "Both BAO and LRG likelihoods should not be used at the same time.");
        if(usePlanck_)
        {
            planckLike_.reset(new PlanckLikelihood(true, true, true, useHighP, planckLikeLite, false, false, false, 100, false));
            if(!planckLikeLite)
                planckLike_->setSZPrior(true);
        }
        if(useBAO_)
        {
            likes_.push_back(new BAOLikelihood(cmb_, false));
        }
        if(useWiggleZ_)
        {
            likes_.push_back(new WiggleZLikelihood(datapath, cmb_, 'a', false)); 
            likes_.push_back(new WiggleZLikelihood(datapath, cmb_, 'b', false)); 
            likes_.push_back(new WiggleZLikelihood(datapath, cmb_, 'c', false)); 
            likes_.push_back(new WiggleZLikelihood(datapath, cmb_, 'd', false)); 
        }
        if(useSZ_)
        {
            likes_.push_back(new PlanckSZLikelihood(cmb_, false));
        }
        if(useH0_)
        {
            likes_.push_back(new H0Likelihood);
        }
    }

    ~CombinedLikelihood()
    {
    }

    void setCosmoParams(const CosmologicalParams& params)
    {
        params_ = &params;
        cmb_.initialize(params, true, true, true, true, 1.0);
        if(usePlanck_)
        {
            cmb_.getLensedCl(&clTT_, &clEE_, &clTE_, &clBB_);
            cmb_.getCl(NULL, NULL, NULL, &clPP_, NULL, NULL);
            planckLike_->setCls(&clTT_, &clEE_, &clTE_, &clBB_, &clPP_);
        }
        for(auto like : likes_)
            like->setCosmoParams(params);
    }

    double likelihood()
    {
        double lnLike = 0; // This is -2*ln(likelihood)
        if(usePlanck_)
            lnLike = lnLike + planckLike_->likelihood();
        for(auto like : likes_)
            lnLike += like->likelihood();
        return lnLike;
    }

    void setModelCosmoParams(CosmologicalParams *params)
    {
        // Sets modelParams_ to initial cosmological parameters.
        modelParams_ = params;
        // Sets vModel_ to a vector containing the parameters in modelParams_.
        modelParams_->getAllParameters(vModel_);

        // set all to 0
        for(auto it = vModel_.begin(); it != vModel_.end(); ++it)
            *it = 0;
    }

    double calculate(double* params, int nPar)
    {
        // Check to see that modelParams_ is set. This is set by calling setModelCosmoParams().
        check(modelParams_, "model params must be set before calling this function");
        // Check to see that vModel_ is not empty
        check(!vModel_.empty(), "");
        const int nModel = vModel_.size();
        
        // Check that number of parameters passed in is same as number of parameters in nModel
        int extraPar = 0;
        if(usePlanck_)
        {
            extraPar = 1;
            if(!planckLikeLite_)
            {
                extraPar += 15;
                if(useHighP_)
                    extraPar += 18;
            }
        }
        check(nPar == nModel + extraPar, "wrong number of model params: expected " << nModel + extraPar << ", provided " << nPar);
    
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

        if(usePlanck_)
        {
            planckLike_->setAPlanck(params[nModel]);
            if(!planckLikeLite_)
            {
                int extraSize = 15 + (useHighP_ ? 17 : 0);
                check(planckExtraParams_.size() == extraSize, "");
                for(int i = 0; i < extraSize; ++i)
                    planckExtraParams_[i] = params[nModel + 1 + i];
                planckLike_->setHighExtraParams(planckExtraParams_);
                if(useHighP_)
                    planckLike_->setAPol(params[nPar - 1]);
            }
        }
    
        return likelihood();
    }

private:
    CMB cmb_;
    const CosmologicalParams* params_; // Cosmological parameters for initialization
    

    std::vector<double> clTT_, clEE_, clTE_, clBB_, clPP_;

    CosmologicalParams* modelParams_; // Cosmological parameters in test models
    std::vector<double> vModel_; // modelParams_ as a vector

    std::vector<double> planckExtraParams_;

    // Specifies which likelihoods to include
    bool usePlanck_, useBAO_, useWiggleZ_, useSZ_, useH0_;
    bool useHighP_, planckLikeLite_;

    // Likelihood objects
    std::unique_ptr<PlanckLikelihood> planckLike_;
    std::vector<Math::CosmoLikelihood*> likes_;
};

#endif

