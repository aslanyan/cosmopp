#ifndef COSMO_PP_PLANCK_LIKE_FAST_HPP
#define COSMO_PP_PLANCK_LIKE_FAST_HPP

#include <vector>

#include <planck_like.hpp>
#include <cmb.hpp>
#include <learn_as_you_go.hpp>

class PlanckLikeFast : public Math::LikelihoodFunction
{
public:
    PlanckLikeFast(CosmologicalParams* params, bool useCommander = true, bool useCamspec = true, bool useLensing = true, bool usePolarization = false, bool useActSpt = false, bool includeTensors = false, double kPerDecade = 100, double precision = 0.2, unsigned long minCount = 10000);
    ~PlanckLikeFast();

    double calculate(double* params, int nPar) { return doCalculation(params, nPar, false); }
    double calculateExact(double* params, int nPar) { return doCalculation(params, nPar, true); }

    void setPrecision(double p);

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
};

#endif

