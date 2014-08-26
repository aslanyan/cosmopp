#ifndef COSMO_PP_MATTER_LIKELIHOOD_HPP
#define COSMO_PP_MATTER_LIKELIHOOD_HPP

#include <vector>

#include <function.hpp>
#include <cosmological_params.hpp>

#include <gmd.h>
#include <lavd.h>

class MatterLikelihood
{
public:
    MatterLikelihood(const char* pkFileName, const char* covFileName);

    double calculate(const Math::RealFunction& matterPk, const CosmologicalParams& params) const;

    double calculateLin(const Math::RealFunction& matterPk, const CosmologicalParams& params, double b2) const;

    void useScaling(const CosmologicalParams& paramsFid, double z);

    double DV(const CosmologicalParams& params, double z) const;

private:
    void readPk(const char* fileName);
    void readLogCov(const char* fileName);
    

private:
    LaGenMatDouble cInv_;
    std::vector<double> kh_, data_;
    bool scale_;
    double dvFid_;
    double z_;
};

#endif

