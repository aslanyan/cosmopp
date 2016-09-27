#pragma once
#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <string>
#include <fstream>

#include <macros.hpp>
#include <planck_like.hpp>
#include <wmap9_like.hpp>
#include <bao_like.hpp>
#include <numerics.hpp>
#include <wigglez_like.hpp>
#include <plancksz_like.hpp>
#include <h0_like.hpp>

class CombinedLikelihood : public Math::LikelihoodFunction
{
public:
    CombinedLikelihood(bool usePlanck, bool useWMAP, bool useBAO, bool useWiggleZ, bool useSZ = false, bool useH0 = false, bool useHighP = false, bool planckLikeLite = true) : usePlanck_(usePlanck), useWMAP_(useWMAP), useBAO_(useBAO), useLRG_(useLRG), useWiggleZ_(useWiggleZ), useSZ_(useSZ), useH0_(useH0), useHighP_(useHighP), planckLikeLite_(planckLikeLite), planckExtraParams_(planckLikeLite ? 0 : 15 + (useHighP ? 17 : 0))
    {
        cosmo_.preInitialize(3500, false, true, false, 0, 100, 1e-6, 1.0);
        nLikes_ = 0;
        //check(!(usePlanck_ && useWMAP_), "Both Planck and WMAP likelihoods should not be used at the same time.");
        //check(!(useBAO_ && useLRG_), "Both BAO and LRG likelihoods should not be used at the same time.");
        if(usePlanck_)
        {
            planckLike_ = new PlanckLikelihood(true, true, true, useHighP, planckLikeLite, false, false, false, 100, false);
            if(!planckLikeLite)
                planckLike_->setSZPrior(true);
        }
        if(useWMAP_)
            wmapLike_ = new WMAP9Likelihood(true, true, true, true, true, true);
        if(useBAO_)
        {
            likes_.push_back(new BAOLikelihood(cosmo_, false, useLRG_));
            ++nLikes_;
        }
        if(useLRG_)
        {
            likes_.push_back(new LRGDR7Likelihood(datapath, cosmo_, false));
            ++nLikes_;
        }
        if(useWiggleZ_)
        {
            likes_.push_back(new WiggleZLikelihood(datapath, cosmo_, 'a', false)); 
            likes_.push_back(new WiggleZLikelihood(datapath, cosmo_, 'b', false)); 
            likes_.push_back(new WiggleZLikelihood(datapath, cosmo_, 'c', false)); 
            likes_.push_back(new WiggleZLikelihood(datapath, cosmo_, 'd', false)); 
            nLikes_ += 4;
        }
        if(useSZ_)
        {
            likes_.push_back(new PlanckSZLikelihood(cosmo_, false));
            ++nLikes_;
        }
        if(useH0_)
        {
            likes_.push_back(new H0Likelihood);
            ++nLikes_;
        }
    }

    ~CombinedLikelihood()
    {
        if(usePlanck_)
            delete planckLike_;
        if(useWMAP_)
            delete wmapLike_;
    }

    void setCosmoParams(const CosmologicalParams& params)
    {
        params_ = &params;
        cosmo_.initialize(params, true, true, true, true, 1.0);
        if(usePlanck_)
        {
            cosmo_.getLensedCl(&clTT_, &clEE_, &clTE_, &clBB_);
            cosmo_.getCl(NULL, NULL, NULL, &clPP_, NULL, NULL);
            planckLike_->setCls(&clTT_, &clEE_, &clTE_, &clBB_, &clPP_);
        }
        if(useWMAP_)
        {
            wmapLike_->setCosmoParams(params);
            wmapLike_->calculateCls();
        }
        for(int i = 0; i < nLikes_; ++i)
            likes_[i]->setCosmoParams(params);
    }

    double likelihood()
    {
        double lnLike = 0; // This is -2*ln(likelihood)
        if(usePlanck_)
            lnLike = lnLike + planckLike_->likelihood();
        if(useWMAP_)
            lnLike = lnLike + wmapLike_->likelihood();
        for(int i = 0; i < nLikes_; ++i)
            lnLike += likes_[i]->likelihood();
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
    Cosmo cosmo_;
    const CosmologicalParams* params_; // Cosmological parameters for initialization
    

    std::vector<double> clTT_, clEE_, clTE_, clBB_, clPP_;

    CosmologicalParams* modelParams_; // Cosmological parameters in test models
    std::vector<double> vModel_; // modelParams_ as a vector

    std::vector<double> planckExtraParams_;

    // Specifies which likelihoods to include
    bool usePlanck_, useWMAP_, useBAO_, useLRG_, useWiggleZ_, useSZ_, useH0_;
    bool useHighP_, planckLikeLite_;

    // Likelihood objects
    PlanckLikelihood* planckLike_;
    WMAP9Likelihood* wmapLike_;
    std::vector<Math::CosmoLikelihood*> likes_;
    //BAOLikelihood* BAOLike_;

    // Number of likelihood objects excluding Planck and WMAP
    int nLikes_;
};
