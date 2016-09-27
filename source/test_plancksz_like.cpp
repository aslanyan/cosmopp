#include <macros.hpp>
#include <cosmological_params.hpp>
#include <plancksz_like.hpp>
#include <numerics.hpp>

int main(int argc, char *argv[])
{
    try {
        CMB cmb;
        cmb.preInitialize(3500, false, true, false);
        PlanckSZLikelihood planckSZ(cmb);
        
        const double ombh2 = 0.022;
        const double omch2 = 0.11;
        const double h = 0.7;
        const double tau = 0.1;

        const double ns = 0.9619;
        const double as = 2.2154e-9;
        const double pivot = 0.05;
        
        LambdaCDMParams params(ombh2, omch2, h, tau, ns, as, pivot);
        planckSZ.setCosmoParams(params);

        const double like = planckSZ.likelihood();

        output_screen("Likelihood = " << like << std::endl);

        if(!Math::areEqual(like, 2.89, 0.1))
        {
            output_screen("FAIL: expected 2.89, got " << like << std::endl);
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
