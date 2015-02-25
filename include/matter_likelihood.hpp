#ifndef COSMO_PP_MATTER_LIKELIHOOD_HPP
#define COSMO_PP_MATTER_LIKELIHOOD_HPP

#include <vector>

#include <function.hpp>
#include <cosmological_params.hpp>
#include <matrix.hpp>

class MatterLikelihood
{
public:
    MatterLikelihood(const char* pkFileName, const char* covFileName, bool isLog = true, double khMin = 0, double khMax = 100);

    double calculate(const Math::RealFunction& matterPk, const CosmologicalParams& params) const;

    double calculateLin(const Math::RealFunction& matterPk, const CosmologicalParams& params) const;

    void useScaling(const CosmologicalParams& paramsFid, double z);

    double DV(const CosmologicalParams& params, double z) const;

private:
    void readPk(const char* fileName, double khMin = 0, double khMax = 100);
    void readLogCov(const char* fileName);
    void readCov(const char* fileName);

private:
    Math::SymmetricMatrix<double> cInv_;
    std::vector<double> kh_, data_;
    bool scale_;
    double dvFid_;
    double z_;

    int nIgnored_, nTotal_;
};

#endif

