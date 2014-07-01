#include <string>
#include <fstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <mcmc.hpp>
#include <markov_chain.hpp>

// Simple two dimensional Gaussian likelihood function
class ExampleMHLikelihood : public Math::LikelihoodFunction
{
public:
    ExampleMHLikelihood() {}
    ~ExampleMHLikelihood() {}

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

int main(int argc, char *argv[])
{
    try {
        StandardException exc;

        using namespace Math;

        // Create the likelihood
        ExampleMHLikelihood like;
        std::string root = "example_files/metropolis_hastings";

        // Create the Metropolis-Hastings sampler
        MetropolisHastings mh(2, like, root);

        // Assign parameter names and ranges
        const double xMin = -10, xMax = 10, yMin = -10, yMax = 10;
        mh.setParam(0, "x", xMin, xMax, 0, 5, 1, 0.05);
        mh.setParam(1, "y", yMin, yMax, 0, 5, 1, 0.05);

        // Choose the burin and run
        const unsigned long burnin = 1000;
        const int nChains = mh.run(100000, 0, burnin, MetropolisHastings::GELMAN_RUBIN, 0.001);

        // Read the resulting chain(s) without thinning
        const unsigned int thin = 1;
        MarkovChain chain(nChains, root.c_str(), burnin, thin);

        // Get the one dimensional marginalized posterior distributions, Gaussian smoothed with a scale of 0.3
        Posterior1D* px = chain.posterior(0, Posterior1D::GAUSSIAN_SMOOTHING, 0.2);
        Posterior1D* py = chain.posterior(1, Posterior1D::GAUSSIAN_SMOOTHING, 0.2);

        // Get the two dimensional posterior distribution, gaussian smoothed with a scale of 0.25
        Posterior2D* pxy = chain.posterior(0, 1, 0.2, 0.2);

        // Write the distributions into text files
        px->writeIntoFile("example_files/mh_px.txt");
        py->writeIntoFile("example_files/mh_py.txt");
        pxy->writeIntoFile("example_files/mh_pxy.txt", 5000);

        // Write the contour levels for the 2D distribution into a text file. This can be used later to make contour plots
        std::ofstream out("example_files/mh_contour_levels.txt");
        if(!out)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Cannot write into file example_files/mh_contour_levels.txt.";
            exc.set(exceptionStr.str());
            throw exc;
        }

        out << pxy->get1SigmaLevel() << std::endl << pxy->get2SigmaLevel() << std::endl;
        out.close();

        // Delete the posterior distributions
        delete px;
        delete py;
        delete pxy;

    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
