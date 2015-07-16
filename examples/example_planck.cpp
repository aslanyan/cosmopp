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

        // Create a special case of cosmological params
        LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);

        // Create the Planck likelihood
#ifdef COSMO_PLANCK_15
        PlanckLikelihood planck(true, true, true, true, true, true, true);
#else
        PlanckLikelihood planck(true, true, true, true, false, true);
#endif

        // Set the cosmological parameters
        planck.setCosmoParams(params);

#ifdef COSMO_PLANCK_15
        // Calculate and print the low-l likelihood
        const double lowLike = planck.lowLike();
        output_screen("Low-l likelihood = " << lowLike << std::endl);

        // Calculate and print the high-l likelihood
        const double highLike = planck.highLike();
        output_screen("High-l likelihood = " << highLike << std::endl);

        // Calculate and print the lensing likelihood
        const double lensingLike = planck.lensingLike();
        output_screen("Lensing likelihood = " << lensingLike << std::endl);
#else
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
#endif

        // Calculate and print the total likelihood
        output_screen(std::endl << "Total likelihood = " << planck.likelihood() << std::endl);
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}

