#include <fstream>
#include <vector>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <cosmological_params.hpp>
#include <cmb.hpp>

int main(int argc, char *argv[])
{
    try {
        StandardException exc;

        // Choose values of the cosmological parameters
        const double h = 0.6704;
        const double omBH2 = 0.022032;
        const double omCH2 = 0.12038;
        const double tau = 0.0925;
        const double ns = 0.9619;
        const double as = 2.2154e-9;
        const double pivot = 0.05;

        const double r = 1e-10;
        const double nt = 0;

        const double nEff = 3.046; 
        const int nMassive = 1;
        const double sumMNu = 0.5;

        // Create an instance of CosmologicalParams. Can use other examples below (just uncomment the appropriate line and comment this line)
        LCDMWithTensorParams params(omBH2, omCH2, h, tau, ns, as, pivot, r, nt, pivot); 
        //LCDMWithDegenerateNeutrinosParams params(omBH2, omCH2, h, tau, ns, as, pivot, nEff, nMassive, sumMNu);
        //LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);

        // Choose lMax and create vectors where the resulting cl values will be written
        int lMax = 3000;
        std::vector<double> clTT, clEE, clTE, clPP, clTP, clEP, clBB, clLensedTT, clLensedEE, clLensedTE, clLensedBB;

        // Create a CMB instance.
        CMB cmb;

        output_screen("Pre-initializing CLASS..." << std::endl);
        // pre-initialize cmb
        cmb.preInitialize(lMax, false, false, true, lMax);
        output_screen("OK" << std::endl);

        output_screen("Initializing CLASS..." << std::endl);
        // initialize cmb
        cmb.initialize(params, true, true, true);
        output_screen("OK" << std::endl);

        // get different cl values
        cmb.getCl(&clTT, &clEE, &clTE, &clPP, &clTP, &clEP, &clBB);
        // get different lensed cl values
        cmb.getLensedCl(&clLensedTT, &clLensedEE, &clLensedTE, &clLensedBB);

        // write results into output files
        std::ofstream out("example_files/cl.txt");
        std::ofstream outLensed("example_files/cl_lensed.txt");

        for(int l = 0; l <= lMax; ++l)
        {
            out << l  << ' ' << clTT[l] << ' ' << clEE[l] << ' ' << clTE[l] << ' ' << clPP[l] << ' ' << clTP[l] << ' ' << clEP[l] << ' ' << clBB[l] << std::endl;
        }

        for(int l = 0; l < clLensedTT.size(); ++l)
            outLensed << l << ' ' << clLensedTT[l]  << ' ' << clLensedEE[l] << ' ' << clLensedTE[l] << ' ' << clLensedBB[l] << std::endl;

        outLensed.close();
        out.close();
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}

