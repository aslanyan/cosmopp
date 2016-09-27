#include <macros.hpp>
#include <bao_like.hpp>
#include <cosmological_params.hpp>
#include <numerics.hpp>

int main(int argc, char *argv[])
{
    try {
        CMB cmb;
        cmb.preInitialize(3500);
        BAOLikelihood bao(cmb);
        
        const double ombh2 = 0.022;
        const double omch2 = 0.11;
        const double h = 0.73;
        const double tau = 0.1;

        const double ns = 0.9619;
        const double as = 2.2154e-9;
        const double pivot = 0.05;
        
        LambdaCDMParams params(ombh2, omch2, h, tau, ns, as, pivot);
        bao.setCosmoParams(params);

        output_screen("  6df likelihood = " << bao.SixDF_Likelihood() << std::endl);
        output_screen(" lowZ likelihood = " << bao.LowZ_DR10_11_Likelihood() << std::endl);
        output_screen("cmass likelihood = " << bao.CMASS_DR10_11_Likelihood() << std::endl);
        output_screen("  mgs likelihood = " << bao.MGS_DR7_Likelihood() << std::endl);

        const double like = bao.likelihood();
        if(!Math::areEqual(like, 47.5, 0.01))
        {
            output_screen("FAIL: expected 47.5, got " << like << std::endl);
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
