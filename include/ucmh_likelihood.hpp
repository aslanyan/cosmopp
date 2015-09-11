#ifndef COSMO_PP_UCMH_LIKELIHOOD_HPP
#define COSMO_PP_UCMH_LIKELIHOOD_HPP

#include <map>
#include <vector>

#include <function.hpp>
#include <table_function.hpp>

class UCMHLikelihood
{
public:
    UCMHLikelihood(const char* fileName, bool lateKineticDecoupling = false);
    ~UCMHLikelihood() {}

    double calculate(const Math::RealFunction& ps) const;

private:
    void extractCLValues(const std::string& s, std::vector<double>& cl) const;
    void createClToLike();

    void readIntegralFactor();

    double calculateSigma(double k, const Math::RealFunction *ps = NULL) const;

private:
    typedef Math::TableFunction<double, double> LikelihoodType;

    std::map<double, LikelihoodType> likes_;
    Math::TableFunction<double, double> clToLike_;

    const double pkMax_, likeMax_;

    Math::TableFunction<double, double> integralFactor_;
    double sigma0_;
};

#endif

