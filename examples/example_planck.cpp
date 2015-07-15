#include <fstream>
#include <vector>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>

#include <planck_like.hpp>

int main(int argc, char *argv[])
{
    try {
        // Choose values of the cosmological parameters
        const double omBH2 = 0.022032;
        const double omCH2 = 0.12038;
        const double h = 0.6704;
        const double tau = 0.0925;
        const double ns = 0.9619;
        const double as = 2.2154e-9;
        const double pivot = 0.05;

        const double r = 0.2;
        const double nt = 0;

        const double nEff = 3.046; 
        const int nMassive = 1;
        const double sumMNu = 0.0;

#ifdef COSMO_PLANCK_15
        output_screen("This is not implemented for Planck 15!" << std::endl);
        return 1;
#else

        // Create a special case of cosmological params
        LCDMWithTensorAndDegenerateNeutrinosParams params(omBH2, omCH2, h, tau, ns, as, pivot, r, nt, pivot, nEff, nMassive, sumMNu);

        // Create the Planck likelihood
        PlanckLikelihood planck(true, true, true, true, false, true);

        // Set the cosmological parameters
        planck.setCosmoParams(params);

        // Set the foreground parameters for Camspec
        planck.setCamspecExtraParams(153, 54.9, 55.8, 4, 55.5, 4, 0.91, 0.63, 0.6, 1, 1, 0.1, 1, 0.3);

        // Calculate and print the commander likelihood
        const double comLike = planck.commanderLike();
        output_screen("Commander likelihood = " << comLike << std::endl);

        // Calculate and print the Camspec likelihood
        const double camspecLike = planck.camspecLike();
        output_screen("Camspec likelihood = " << camspecLike << std::endl);

        // Calculate and print the WMAP polarization likelihood
        const double polLike = planck.polLike();
        output_screen("Pol likelihood = " << polLike << std::endl);

        // Calculate and print the lensing likelihood
        const double lensingLike = planck.lensingLike();
        output_screen("Lensing likelihood = " << lensingLike << std::endl);

        // Calculate and print the total likelihood
        output_screen(std::endl << "Total likelihood = " << planck.likelihood() << std::endl);
#endif
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}

