#include <fstream>
#include <string>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <likelihood_function.hpp>
#include <mn_scanner.hpp>
#include <planck_like.hpp>

int main(int argc, char *argv[])
{
    try {
        StandardException exc;

#ifdef COSMO_PLANCK_15
        // Create the likelihood
        PlanckLikelihood like(true, true, true, false, true, false, false);

        // Create the Multinest scanner with 7 parameters and 300 live points
        MnScanner scanner(7, like, 300, std::string("example_files/planck_multinest_"));
#else
        // Create the likelihood
        PlanckLikelihood like(true, true, true, true);

        // Create the Multinest scanner with 20 parameters and 300 live points
        MnScanner scanner(20, like, 300, std::string("example_files/planck_multinest_"));
#endif

        // Set the parameter names and ranges
        scanner.setParam(0, std::string("ombh2"), 0.02, 0.025);
        scanner.setParam(1, std::string("omch2"), 0.1, 0.2);
        scanner.setParam(2, std::string("h"), 0.55, 0.85);
        scanner.setParam(3, std::string("tau"), 0.01, 0.30);
        scanner.setParam(4, std::string("ns"), 0.9, 1.1);
        scanner.setParam(5, std::string("as"), 2.7, 3.5);

#ifdef COSMO_PLANCK_15
        scanner.setParamGauss(6, "A_planck", 1.0, 0.0025);
#else
        scanner.setParam(6, "A_ps_100", 0, 360);
        scanner.setParam(7, "A_ps_143", 0, 270);
        scanner.setParam(8, "A_ps_217", 0, 450);
        scanner.setParam(9, "A_cib_143", 0, 20);
        scanner.setParam(10, "A_cib_217", 0, 80);
        scanner.setParam(11, "A_sz", 0, 10);
        scanner.setParam(12, "r_ps", 0.0, 1.0);
        scanner.setParam(13, "r_cib", 0.0, 1.0);
        scanner.setParam(14, "n_Dl_cib", -2, 2);
        scanner.setParam(15, "cal_100", 0.98, 1.02);
        scanner.setParam(16, "cal_127", 0.95, 1.05);
        scanner.setParam(17, "xi_sz_cib", 0, 1);
        scanner.setParam(18, "A_ksz", 0, 10);
        scanner.setParam(19, "Bm_1_1", -20, 20);
#endif

        // Set the model for planck likelihood by specifying an example parameter set.
        const double pivot = 0.05;
        LambdaCDMParams par(0.022, 0.12, 0.7, 0.1, 1.0, std::exp(3.0) / 1e10, pivot);
        like.setModelCosmoParams(&par);

        // Run the scanner. The Results will be output in corresponding files at the end
        scanner.run();
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}

