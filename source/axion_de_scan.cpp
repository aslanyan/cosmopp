#include <cosmo_mpi.hpp>

#include <cmath>

#include <cosmological_params.hpp>
#include <math_constants.hpp>
#include <mcmc.hpp>
#include <macros.hpp>
#include <planck_like.hpp>
#include <markov_chain.hpp>
#include <timer.hpp>

namespace
{

class AxionDEParamsPhenom : public LambdaCDMParams
{
public:
    AxionDEParamsPhenom(double omBH2, double omLambdaH2, double h, double tau, double ns, double as, double pivot, double nf) : LambdaCDMParams(omBH2, h * h - omBH2 - omLambdaH2, h, tau, ns, as, pivot), nf_(nf), omLambdaH2_(omLambdaH2)
    {
    }

    virtual std::string name() const { return "AxionDEParamsPhenom"; }

    // the parameters are 
    virtual void getAllParameters(std::vector<double>& v) const
    {
        v.resize(7);
        v[0] = getOmBH2();
        v[1] = omLambdaH2_;
        v[2] = getH();
        v[3] = getTau();
        v[4] = getNs();
        v[5] = std::log(getAs() * 1e10);
        v[6] = nf_;
    }

    virtual bool setAllParameters(const std::vector<double>& v, double *badLike = NULL)
    {
        check(v.size() >= 7, "");
        omBH2_ = v[0];
        omLambdaH2_ = v[1];
        h_ = v[2];
        tau_ = v[3];
        ps_.setNs(v[4]);
        ps_.setAs(std::exp(v[5]) / 1e10);
        nf_ = v[6];

        omCH2_ = h_ * h_ - omBH2_ - omLambdaH2_;

        if(badLike)
            *badLike = 0;

        return true;
    }

private:
    double nf_;
    double omLambdaH2_;
};

class AxionDEPriorPhenom : public Math::PriorFunctionBase
{
public:
    AxionDEPriorPhenom(double nfMean) : p_(1.0 / nfMean)
    {
        check(nfMean >= 1, "");
        check(p_ > 0 && p_ < 1, "");
    }

    virtual double calculate(double* params, int nPar)
    {
        check(nPar == 8, "");
        double res = 1;

        const double omBH2 = params[0];
        const double omLambdaH2 = params[1];
        const double h = params[2];
        const double tau = params[3];
        const double ns = params[4];
        const double as = params[5];
        const double nf = params[6];
        const double Aplanck = params[7];

        const double omBH2Min = 0.005, omBH2Max = 0.1;
        const double hMin = 0.2, hMax = 1.0;
        const double tauMin = 0.01, tauMax = 0.8;
        const double nsMin = 0.9, nsMax = 1.1;
        const double asMin = 2.7, asMax = 4.0;
        const double AplanckMean = 1.0, AplanckSigma = 0.0025;

        if(omBH2 < omBH2Min || omBH2 > omBH2Max)
            return 0;
        res *= 1.0 / (omBH2Max - omBH2Min);
        if(h < hMin || h > hMax)
            return 0;
        res *= 1.0 / (hMax - hMin);
        if(tau < tauMin || tau > tauMax)
            return 0;
        res *= 1.0 / (tauMax - tauMin);
        if(ns < nsMin || ns > nsMax)
            return 0;
        res *= 1.0 / (nsMax - nsMin);
        if(as < asMin || as > asMax)
            return 0;
        res *= 1.0 / (asMax - asMin);

        const double deltaAplanck = Aplanck - AplanckMean;
        res *= (std::exp(-deltaAplanck * deltaAplanck / (2 * AplanckSigma * AplanckSigma)) / (AplanckSigma * std::sqrt(2 * Math::pi)));

        // geometric distrib small p approximation
        if(nf < 1)
            return 0;
        res *= (p_ * std::pow(1.0 - p_, nf - 1));

        const double rhoLambda = omLambdaH2 * 8.094e-11; // ev^4

        const double muM = 1.277e-68; // ev^2
        const double sigmaM = 7.889e-68; // ev^2
        const double muPhi = 1.976e54; // ev^2
        const double sigmaPhi = 1.768e54; // ev^2
        const double muLambda = nf * muM * muPhi / 2.0;
        const double sigmaLambda2 = (sigmaM * sigmaM * sigmaPhi * sigmaPhi + sigmaM * sigmaM * muPhi * muPhi + sigmaPhi * sigmaPhi * muM * muM) * nf / 2.0;

        // rhoLambdaPrior
        const double deltaRhoLambda = rhoLambda - muLambda;
        res *= (std::exp(-deltaRhoLambda * deltaRhoLambda / (2 * sigmaLambda2)) / std::sqrt(2 * Math::pi * sigmaLambda2));

        return res;
    }

private:
    const double p_;
};

} // namespace

int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        if(argc < 2)
        {
            std::string exceptionStr = "N_f mean must be specified.";
            exc.set(exceptionStr);
            throw exc;
        }

        double nfMean;

        std::stringstream inStr;
        inStr << argv[1];
        inStr >> nfMean;

        if(nfMean < 1 || nfMean > 1e9)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Invalid value of " << nfMean << " for N_f mean.";
            exc.set(exceptionStr.str());
            throw exc;
        }

        using namespace Math;

        std::string root = "axion_de";
#ifdef COSMO_PLANCK_15
        PlanckLikelihood planckLike(true, true, true, true, true, false, false, false, 5);
        MetropolisHastings mh(8, planckLike, root);

        mh.setParam(0, "ombh2", 0.005, 0.1, 0.022, 0.0003, 0.00005);
        mh.setParam(1, "omlambdah2", 0.01, 0.9, 0.34, 0.003, 0.0005);
        mh.setParam(2, "h", 0.2, 1.0, 0.7, 0.02, 0.002);
        mh.setParam(3, "tau", 0.01, 0.8, 0.1, 0.01, 0.002);
        mh.setParam(4, "ns", 0.9, 1.1, 1.0, 0.01, 0.002);
        mh.setParam(5, "As", 2.7, 4.0, 3.0, 0.1, 0.002);
        mh.setParam(6, "nf", 1.0, nfMean * nfMean, nfMean, nfMean / 5, nfMean / 25);
        mh.setParamGauss(7, "A_planck", 1.0, 0.0025, 1.0, 0.001, 0.0002);

        AxionDEPriorPhenom prior(nfMean);
        mh.useExternalPrior(&prior);

        const double pivot = 0.05;
        AxionDEParamsPhenom par(0.022, 0.34, 0.7, 0.1, 1.0, std::exp(3.0) / 1e10, pivot, nfMean);
        planckLike.setModelCosmoParams(&par);

        Timer timer("MCMC PLANCK AXION DE");

        const unsigned long burnin = 1000;
        timer.start();
        const int nChains = mh.run(50000, 1, burnin, MetropolisHastings::GELMAN_RUBIN, 0.01, true);
        const unsigned long time = timer.end();
        output_screen("MCMC Planck Axion DE took " << time / 1000000 << " seconds." << std::endl);

        if(CosmoMPI::create().isMaster())
            return 0;
        
        const unsigned int thin = 2;
        MarkovChain chain(nChains, root.c_str(), burnin, thin);

        const int nPoints = 1000;
        const int nPar = 8;

        std::ofstream outParamLimits("axion_de_param_limits.txt");
        for(int i = 0; i < nPar; ++i)
        {
            const std::string& paramName = mh.getParamName(i);
            std::stringstream fileName;
            fileName << "axion_de_" << paramName << ".txt";
            Posterior1D* p = chain.posterior(i, Posterior1D::GAUSSIAN_SMOOTHING);

            p->writeIntoFile(fileName.str().c_str(), nPoints);

            const double median = p->median();
            double lower, upper;
            p->get1SigmaTwoSided(lower, upper);
            const double sigma = (upper - lower) / 2.0;

            outParamLimits << paramName << " = " << median << "+-" << sigma << std::endl;
        }
        outParamLimits.close();
#else
        output_screen("THIS IS ONLY IMPLEMENTED FOR PLANCK LIKELIHOOD 15!" << std::endl);
#endif
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
