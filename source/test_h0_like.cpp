#include <macros.hpp>
#include <cosmological_params.hpp>
#include <h0_like.hpp>
#include <numerics.hpp>

int main(int argc, char *argv[])
{
    try {
        H0Likelihood H0(73.03, 1.79);
        
        const double ombh2 = 0.022;
        const double omch2 = 0.11;
        const double h = 0.67;
        const double tau = 0.1;

        const double ns = 0.9619;
        const double as = 2.2154e-9;
        const double pivot = 0.05;
        
        LambdaCDMParams params(ombh2, omch2, h, tau, ns, as, pivot);
        H0.setCosmoParams(params);

        const double like = H0.likelihood();

        output_screen("Likelihood = " << like << std::endl);

        if(!Math::areEqual(like, 11.3482, 0.01))
        {
            output_screen("FAIL: expected 11.3482, got " << like << std::endl);
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
