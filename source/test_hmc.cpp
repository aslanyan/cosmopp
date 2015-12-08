#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <memory>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <hmc.hpp>
#include <markov_chain.hpp>

namespace
{

class UncorrelatedMultiDimLike : public Math::LikelihoodWithDerivs
{
public:
    UncorrelatedMultiDimLike(const std::vector<double>& mean, const std::vector<double>& sigma) : mean_(mean), sigma_(sigma)
    {
        check(!mean_.empty(), "");
        check(mean_.size() == sigma_.size(), "");
    }

    ~UncorrelatedMultiDimLike() {}

    virtual double calculate(double* params, int nParams)
    {
        check(nParams == mean_.size(), "");
        check(mean_.size() == sigma_.size(), "");

        double res = 0;
        for(int i = 0; i < mean_.size(); ++i)
        {
            const double delta = params[i] - mean_[i];
            res += delta * delta / (sigma_[i] * sigma_[i]);
        }

        return res;
    }

    virtual double calculateDeriv(double *params, int nParams, int i)
    {
        check(nParams == mean_.size(), "");
        check(mean_.size() == sigma_.size(), "");
        check(i >= 0 && i < nParams, "");

        return 2 * (params[i] - mean_[i]) / (sigma_[i] * sigma_[i]);
    }

private:
    std::vector<double> mean_, sigma_;
};

// Simple two dimensional Gaussian likelihood function
class ExampleHMCLikelihood : public Math::LikelihoodFunction
{
public:
    ExampleHMCLikelihood() {}
    ~ExampleHMCLikelihood() {}

    virtual double calculate(double* params, int nParams)
    {
        check(nParams == 2, "");
        const double x = params[0], y = params[1];

        const double x1 = (x + y) / 2, y1 = (x - y) / 2;
        const double x0 = 0, y0 = 0;
        const double sigmaX = 1, sigmaY = 2;
        const double deltaX = x1 - x0;
        const double deltaY = y1 - y0;

        return deltaX * deltaX / (sigmaX * sigmaX) + deltaY * deltaY / (sigmaY * sigmaY);
    }
};

} // namespace

int main(int argc, char *argv[])
{
    try {
        StandardException exc;

        using namespace Math;

        const int n = 5;
        std::vector<double> mean(n, 0);
        std::vector<double> sigma(n, 5);

        std::vector<double> mass(n, 2);

        UncorrelatedMultiDimLike like(mean, sigma);
        std::string root = "test_files/hmc_test";

        HMC hmc(n, like, root, 1.0, 10);

        for(int i = 0; i < n; ++i)
        {
            std::stringstream paramName;
            paramName << "x_" << i;
            hmc.setParam(i, paramName.str(), mass[i], 5);
        }

        hmc.run(10000);

        const unsigned int thin = 1;
        const unsigned long burnin = 10;
        MarkovChain chain(1, root.c_str(), burnin, thin);
        for(int i = 0; i < n; ++i)
        {
            std::unique_ptr<Posterior1D> p(chain.posterior(i, Posterior1D::GAUSSIAN_SMOOTHING));
            std::stringstream fileName;
            fileName << "test_files/hmc_" << hmc.getParamName(i) << ".txt";
            p->writeIntoFile(fileName.str().c_str());

            const double median = p->median();
            double l, u;
            p->get1SigmaTwoSided(l, u);
            output_screen("Param " << i << ":\t" << median << " + " << u - median << " - " << median - l << std::endl);
        }
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}

