#ifndef COSMO_PP_H0_LIKE_HPP
#define COSMO_PP_H0_LIKE_HPP

#include <string>
#include <fstream>
#include <vector>

#include <macros.hpp>
#include <likelihood_function.hpp>

/// H0 Likelihood (written by Nicolas Canac)
class H0Likelihood : public Math::CosmoLikelihood
{
public:
    /// Constructor.
    /// \param mean The mean value of H_0. The default value is from http://arxiv.org/abs/1604.01424.
    /// \param stdev The standard deviation of H_0. The default value is from http://arxiv.org/abs/1604.01424.
    H0Likelihood(double mean = 73.03, double stdev = 1.79) : mean_(mean), stdev_(stdev) {}

    /// Destructor.
    ~H0Likelihood() {}

    /// Set the cosmological parameters. This sets the parameters for which the bao likelihoods will be subsequently calculated.
    /// \param params The cosmological parameters object.
    void setCosmoParams(const CosmologicalParams& params)
    {
        params_ = &params;
    }

    /// Likelihood.
    /// \return -2ln(likelihood).
    double likelihood()
    {
        const double H = params_->getH() * 100;
        const double LnLike = pow(H - mean_, 2) / pow(stdev_, 2);
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
    const double mean_, stdev_;
    const CosmologicalParams* params_;
    
    CosmologicalParams* modelParams_;
    std::vector<double> vModel_; // modelParams_ as a vector
};

#endif

