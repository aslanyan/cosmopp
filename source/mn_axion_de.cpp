#include <cosmo_mpi.hpp>

#include <cmath>

#include <cosmological_params.hpp>
#include <math_constants.hpp>
#include <mn_scanner.hpp>
#include <macros.hpp>
#include <planck_like.hpp>
#include <markov_chain.hpp>
#include <timer.hpp>

namespace
{

class AxionDELikelihood : public Math::LikelihoodFunction
{
public:
    AxionDELikelihood() : planck_(true, true, true, true, true, false, false, false, 5), par_(0.022, 0.12, 0.7, 0.1, 1.0, std::exp(3.0) / 1e10, 0.05)
    {
        planck_.setModelCosmoParams(&par_);
    }
    ~AxionDELikelihood() {}

    virtual double calculate(double* params, int nParams)
    {
        check(nParams == 8, "");

        const double omBH2 = params[0];
        const double omLambdaH2 = params[1];
        const double h = params[2];
        const double omCH2 = h * h - omBH2 - omLambdaH2;
        check(omCH2 / (h * h) < 1, "");

        // bad likelihood
        if(omCH2 < 0)
            return (1 - omCH2) * 1e10;

        const double tau = params[3];
        const double ns = params[4];
        const double as = params[5];
        const double nf = params[6];
        const double A_planck = params[7];

        params_[0] = omBH2;
        params_[1] = omCH2;
        params_[2] = h;
        params_[3] = tau;
        params_[4] = ns;
        params_[5] = as;
        params_[6] = A_planck;

        double res = planck_.calculate(params_, 7);

        const double muM = 1.277e-68; // ev^2
        const double sigmaM = 7.889e-68; // ev^2
        const double muPhi = 1.976e54; // ev^2
        const double sigmaPhi = 1.768e54; // ev^2
        const double muLambda = nf * muM * muPhi / 2.0;
        const double sigmaLambda2 = (sigmaM * sigmaM * sigmaPhi * sigmaPhi + sigmaM * sigmaM * muPhi * muPhi + sigmaPhi * sigmaPhi * muM * muM) * nf / 2.0;

        const double rhoLambda = omLambdaH2 * 8.094e-11; // ev^4
        const double deltaRhoLambda = rhoLambda - muLambda;

        res += (std::log(sigmaLambda2) + deltaRhoLambda * deltaRhoLambda / sigmaLambda2);
        
        return res;
    }

private:
    LambdaCDMParams par_;
    PlanckLikelihood planck_;
    double params_[7];
};

class NfPrior : public Math::RealFunction
{
public:
    NfPrior(double nfMean) : p_(1.0 / nfMean)
    {
        check(nfMean >= 1, "");
    }

    virtual double evaluate(double x) const
    {
        check(x >= 1, "");
        return p_ * std::pow(1.0 - p_, x - 1);
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

        const NfPrior nfPrior(nfMean);

        using namespace Math;

        std::string root = "mn_axion_de";
        AxionDELikelihood like;
        MnScanner mn(8, like, 300, root);

        mn.setParam(0, "ombh2", 0.02, 0.025);
        mn.setParam(1, "omlambdah2", 0.2, 0.5);
        mn.setParam(2, "h", 0.55, 0.85);
        mn.setParam(3, "tau", 0.02, 0.20);
        mn.setParam(4, "ns", 0.9, 1.1);
        mn.setParam(5, "As", 2.7, 3.5);
        mn.setParamGeneral(6, "nf", 1.0, 100 * nfMean, nfPrior);
        mn.setParamGauss(7, "A_planck", 1.0, 0.0025);

        Timer timer("MN PLANCK AXION DE");

        timer.start();
        mn.run(true);
        const unsigned long time = timer.end();
        output_screen("MN Planck Axion DE took " << time / 1000000 << " seconds." << std::endl);

        if(CosmoMPI::create().isMaster())
            return 0;
        
        MarkovChain chain("mn_axion_de.txt");

        const int nPoints = 1000;
        const int nPar = 8;

        std::ofstream outParamLimits("mn_axion_de_param_limits.txt");
        for(int i = 0; i < nPar; ++i)
        {
            const std::string& paramName = mn.getParamName(i);
            std::stringstream fileName;
            fileName << "mn_axion_de_" << paramName << ".txt";
            Posterior1D* p = chain.posterior(i, Posterior1D::GAUSSIAN_SMOOTHING);

            p->writeIntoFile(fileName.str().c_str(), nPoints);

            const double median = p->median();
            double lower, upper;
            p->get1SigmaTwoSided(lower, upper);
            const double sigma = (upper - lower) / 2.0;

            outParamLimits << paramName << " = " << median << "+-" << sigma << std::endl;
        }
        outParamLimits.close();
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}

