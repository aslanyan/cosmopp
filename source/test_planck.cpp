#include <fstream>
#include <vector>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <timer.hpp>

#include <planck_like.hpp>

int main(int argc, char *argv[])
{
    try {
        StandardException exc;

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

        const double kCut = 0.0001;

        Timer t1("CL_CALCULATE");


        //LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);
        //LCDMWithTensorParams params(omBH2, omCH2, h, tau, ns, as, pivot, r, nt, pivot); 
        LCDMWithCutoffTensorDegenerateNeutrinosParams params(omBH2, omCH2, h, tau, kCut, ns, as, pivot, r, nt, pivot, nEff, nMassive, sumMNu);

        PlanckLikelihood planck(true, true, true, true, false, true);
        planck.setCosmoParams(params);

        planck.setCamspecExtraParams(153, 54.9, 55.8, 4, 55.5, 4, 0.91, 0.63, 0.6, 1, 1, 0.1, 1, 0.3);
        //planck.setActSptExtraParams(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

        t1.start();
        planck.calculateCls();
        t1.end();

        const double comLike = planck.commanderLike();
        output_screen("Commander likelihood = " << comLike << std::endl);

        const double camspecLike = planck.camspecLike();
        output_screen("Camspec likelihood = " << camspecLike << std::endl);

        const double polLike = planck.polLike();
        output_screen("Pol likelihood = " << polLike << std::endl);

        const double lensingLike = planck.lensingLike();
        output_screen("Lensing likelihood = " << lensingLike << std::endl);

        output_screen(std::endl << "Total likelihood = " << planck.likelihood() << std::endl);
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}

