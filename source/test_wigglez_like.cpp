#include <macros.hpp>
#include <cosmological_params.hpp>
#include <wigglez_like.hpp>
#include <numerics.hpp>

int main(int argc, char *argv[])
{
    try {
        CMB cmb;
        cmb.preInitialize(3500, false, true, false);
        std::string datapath = "./data";

        WiggleZLikelihood likeA(datapath, cmb, 'a');
        WiggleZLikelihood likeB(datapath, cmb, 'b');
        WiggleZLikelihood likeC(datapath, cmb, 'c');
        WiggleZLikelihood likeD(datapath, cmb, 'd');

        const double ombh2 = 0.022;
        const double omch2 = 0.11;
        const double h = 0.73;
        const double tau = 0.1;

        const double ns = 0.9619;
        const double as = 2.2154e-9;
        const double pivot = 0.05;
        
        LCDMWithDegenerateNeutrinosParams params(ombh2, omch2, h, tau, ns, as, pivot, 3.00641, 3, 3 * 0.2);
        likeA.setCosmoParams(params);
        likeB.setCosmoParams(params);
        likeC.setCosmoParams(params);
        likeD.setCosmoParams(params);

        const double like1 = likeA.likelihood();
        const double like2 = likeB.likelihood();
        const double like3 = likeC.likelihood();
        const double like4 = likeD.likelihood();

        output_screen("Likelihood A = " << like1 << std::endl);
        output_screen("Likelihood B = " << like2 << std::endl);
        output_screen("Likelihood C = " << like3 << std::endl);
        output_screen("Likelihood C = " << like4 << std::endl);

        const double totalLike = like1 + like2 + like3 + like4;
        if(!Math::areEqual(totalLike, 480.48, 0.001))
        {
            output_screen("FAIL: expected 480.48, got " << totalLike << std::endl);
            return 1;
        }
        else
        {
            output_screen("SUCCESS!" << std::endl);
        }
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
