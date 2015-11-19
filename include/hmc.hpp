#ifndef COSMO_PP_HMC_HPP
#define COSMO_PP_HMC_HPP

#include <vector>
#include <string>
#include <limits>
#include <ctime>
#include <memory>
#include <fstream>
#include <random>

#include <macros.hpp>
#include <likelihood_function.hpp>

namespace Math
{

class HMC
{
public:
    HMC(int nPar, LikelihoodWithDerivs& like, std::string fileRoot, double tauMax, int nMax, time_t seed = 0);
    ~HMC();

    void setParam(int i, const std::string& name, double mass, double starting);

    /// Get the name of a parameter.
    /// \param i The index of the parameter.
    /// \return The name of the parameter.
    const std::string& getParamName(int i) const { check(i >= 0 && i < n_, "invalid index " << i); return paramNames_[i]; }

    void run(int iters = 100);

private:
    void outputCurrent();
    void generateP();
    void calculateDerivs();
    double calculatePLike();

private:
    int n_;
    Math::LikelihoodWithDerivs& like_;
    std::string fileRoot_;
    double tauMax_;
    int nMax_;

    std::vector<double> mass_;
    std::vector<double> starting_;

    std::vector<std::string> paramNames_;

    std::mt19937 gen_;
    std::uniform_real_distribution<double> uniformDist_;
    std::normal_distribution<> gaussDist_;

    std::vector<double> currentPoint_;
    std::vector<double> x_;
    std::vector<double> derivs_;
    std::vector<double> p_;

    double currentLike_;

    std::ofstream outChain_;
};

} // namespace Math

#endif

