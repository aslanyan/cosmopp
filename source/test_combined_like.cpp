#include <fstream>
#include <vector>
#include <sstream>

#include <macros.hpp>
#include <combined_like.hpp>
#include <numerics.hpp>

int main(int argc, char *argv[])
{
    // Choose values of the cosmological parameters
    const double omBH2 = 0.022;
    const double omCH2 = 0.11;
    const double h = 0.7;
    const double tau = 0.1;
    const double ns = 0.9619;
    const double as = 2.2154e-9;
    const double pivot = 0.05;
    
    // Create cosmological params
    LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);

    // Create likelihoods
    bool usePlanck = true;
    bool useBAO = true;
    bool useWiggleZ = true;
    bool useSZ = false;
    bool useH0 = false;

    std::string datapath = "./data";
    CombinedLikelihood like(datapath, usePlanck, useBAO, useWiggleZ, useSZ, useH0, true, true);
    
    // Set the cosmological parameters
    like.setCosmoParams(params);

    // Calculate likelihoods
    double lnlike = like.likelihood();

    // Output the likelihoods
    output_screen("Combined likelihood = " << lnlike << std::endl);

    const double totalLike = lnlike;

    if(!Math::areEqual(totalLike, 11930.7, 0.0001))
    {
        output_screen("FAIL: expected 11930.7, got " << totalLike << std::endl);
        return 1;
    }
    else
    {
        output_screen("SUCCESS!" << std::endl);
    }

    return 0;
}
