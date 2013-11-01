#include <fstream>
#include <vector>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <phys_constants.hpp>

#include <cosmological_params.hpp>
#include <scale_factor.hpp>
#include <cmb.hpp>

int main(int argc, char *argv[])
{
    try {
        StandardException exc;

        //const double h = 0.6704;
        const double h = 0.611339792294414707818361875979;
        const double omBH2 = 0.0210097173927527372527812588032;
        //const double omBH2 = 0.022032;
        const double omCH2 = 0.135967162847122902125107657412;
        //const double omCH2 = 0.12038;
        const double tau = 0.0366800252574301238750464904115;
        //const double tau = 0.0925;
        const double ns = 0.9619;
        const double as = 2.2154e-9;
        const double pivot = 0.05;

        const double r = 0.1;
        const double nt = 0.0;

        const double nEff = 3.046; 
        const int nMassive = 1;
        const double sumMNu = 0.5;

        std::vector<double> kVals, amplitudes;
        kVals.push_back(1e-6);
        amplitudes.push_back(3e-9);
        kVals.push_back(0.001);
        amplitudes.push_back(1e-9);
        kVals.push_back(3);
        amplitudes.push_back(1e-8);

        //LinearSplineParams params(omBH2, omCH2, h, tau, kVals, amplitudes);
        //LCDMWithDegenerateNeutrinosParams params(omBH2, omCH2, h, tau, ns, as, pivot, nEff, nMassive, sumMNu);
        LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);

        ScaleFactorFunctionClass scaleFactor;
        scaleFactor.initialize(params);

        output_screen("The age of the universe is " << Phys::secToYear(Phys::unitlessToSec(scaleFactor.age())) << " years." << std::endl);
        output_screen("Z_eq = " << params.getOmM() / params.getOmR() - 1 << std::endl);

        int lMax = 3000;
        std::vector<double> clTT, clEE, clTE, clPP, clTP, clEP, clBB, clLensedTT, clLensedEE, clLensedTE, clLensedBB;

        output_screen("Trying out CLASS..." << std::endl);
        CMB cmb;

        output_screen("Pre-initializing CLASS..." << std::endl);
        cmb.preInitialize(lMax, false, true, false, lMax);
        output_screen("OK" << std::endl);

        output_screen("Initializing CLASS..." << std::endl);
        cmb.initialize(params, true, true, true);
        output_screen("OK" << std::endl);

        cmb.getCl(&clTT, &clEE, &clTE, &clPP, &clTP, &clEP, &clBB);
        cmb.getLensedCl(&clLensedTT, &clLensedEE, &clLensedTE, &clLensedBB);

        std::ofstream out("test_cl.txt");
        std::ofstream outLensed("test_cl_lensed.txt");

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

