#ifdef COSMO_MPI
#include <mpi.h>
#endif

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

        // Check if this is the master MPI process
        bool isMaster = true;
#ifdef COSMO_MPI
        int hasMpiInitialized;
        MPI_Initialized(&hasMpiInitialized);
        if(!hasMpiInitialized)
            MPI_Init(NULL, NULL);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(rank != 0)
            isMaster = false;
#endif

        // Create the likelihood
        ExampleMHLikelihood like;
        std::string root = "example_files/metropolis_hastings";

        // Create the Metropolis-Hastings sampler
        MetropolisHastings mh(2, like, root);

        // Assign parameter names and ranges
        const double xMin = -10, xMax = 10, yMin = -10, yMax = 10;
        mh.setParam(0, "x", xMin, xMax, 0, 5, 0.5, 0.05);
        mh.setParam(1, "y", yMin, yMax, 0, 5, 0.5, 0.05);

        // Set the Metropolis-Hastings sampler to vary both parameters together as one block
        std::vector<int> blocks(1, 2);
        mh.specifyParameterBlocks(blocks);

        // Choose the burnin and run
        const unsigned long burnin = 1000;
        const int nChains = mh.run(1000000, 0, burnin, MetropolisHastings::GELMAN_RUBIN, 0.00001);

        // Only the master process will analyze the results
        if(isMaster)
        {
            // Read the resulting chain(s) with thinning
            const unsigned int thin = 1;
            MarkovChain chain(nChains, root.c_str(), burnin, thin);

            // Get the one dimensional marginalized posterior distributions, Gaussian smoothed with a scale of 0.3
            Posterior1D* px = chain.posterior(0, Posterior1D::GAUSSIAN_SMOOTHING);
            Posterior1D* py = chain.posterior(1, Posterior1D::GAUSSIAN_SMOOTHING);

            // Get the two dimensional posterior distribution, gaussian smoothed with a scale of 0.25
            Posterior2D* pxy = chain.posterior(0, 1);

            // Write the distributions into text files

            output_screen("Writing the distributions into text files..." << std::endl);
            px->writeIntoFile("example_files/mh_px.txt");
            py->writeIntoFile("example_files/mh_py.txt");
            pxy->writeIntoFile("example_files/mh_pxy.txt");
            output_screen("OK" << std::endl);

            // Write the contour levels for the 2D distribution into a text file. This can be used later to make contour plots
            std::ofstream out("example_files/mh_contour_levels.txt");
            if(!out)
            {
                std::stringstream exceptionStr;
                exceptionStr << "Cannot write into file example_files/mh_contour_levels.txt.";
                exc.set(exceptionStr.str());
                throw exc;
            }

            out << pxy->get1SigmaLevel() << std::endl;
            //out << pxy->get2SigmaLevel() << std::endl;
            out.close();

            // Delete the posterior distributions
            delete px;
            delete py;
            delete pxy;
        }
        // Finalize MPI if not done so yet
#ifdef COSMO_MPI
        int hasMpiFinalized;
        MPI_Finalized(&hasMpiFinalized);
        if(!hasMpiFinalized)
            MPI_Finalize();
#endif
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
