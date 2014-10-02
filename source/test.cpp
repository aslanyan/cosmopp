#include <cosmo_mpi.hpp>

#include <fstream>
#include <set>
#include <sstream>
#include <iostream>

#include <macros.hpp>
#include <exception_handler.hpp>

#include <test_framework.hpp>
#include <test_unit_conversions.hpp>
#include <test_int_operations.hpp>
#include <test_integral.hpp>
#include <test_conjugate_gradient.hpp>
#include <test_polynomial.hpp>
#include <test_legendre.hpp>
#include <test_spherical_harmonics.hpp>
#include <test_mcmc.hpp>
#include <test_multinest.hpp>
#include <test_mcmc_planck.hpp>
#include <test_multinest_planck.hpp>
#include <test_cmb.hpp>
#include <test_cmb_gibbs.hpp>
#include <test_fit.hpp>
#include <test_planck_like.hpp>
#include <test_wmap9_like.hpp>
#include <test_like_high.hpp>
#include <test_like_low.hpp>
#include <test_wigner_3j.hpp>
#include <test_table_function.hpp>
#include <test_cubic_spline.hpp>
#include <test_three_rotation.hpp>
#include <test_mask_apodizer.hpp>
#include <test_matter_likelihood.hpp>
#include <test_k_nearest_neighbors.hpp>

TestFramework* createTest(const std::string& name)
{
    TestFramework* test = NULL;

    if(name == "unit_conversions")
        test = new TestUnitConversions;
    else if(name == "int_operations")
        test = new TestIntOperations;
    else if(name == "integral")
        test = new TestIntegral;
    else if(name == "conjugate_gradient")
        test = new TestConjugateGradient;
    else if(name == "polynomial")
        test = new TestPolynomial;
    else if(name == "legendre")
        test = new TestLegendre;
    else if(name == "spherical_harmonics")
        test = new TestSphericalHarmonics;
#ifdef COSMO_LAPACKPP
    else if(name == "mcmc_fast")
        test = new TestMCMCFast;
#endif
#ifdef COSMO_MULTINEST
    else if(name == "multinest_fast")
        test = new TestMultinestFast;
#endif
#ifdef COSMO_PLANCK
#ifdef COSMO_CLASS
#ifdef COSMO_LAPACKPP
    else if(name == "mcmc_planck")
        test = new TestMCMCPlanck;
#endif
#ifdef COSMO_MULTINEST
    else if(name == "multinest_planck")
        test = new TestMultinestPlanck;
#endif
#endif
#endif
#ifdef COSMO_CLASS
    else if(name == "cmb")
        test = new TestCMB;
#endif
#ifdef COSMO_HEALPIX
#ifdef COSMO_LAPACKPP
#ifdef COSMO_CLASS
    else if(name == "cmb_gibbs")
        test = new TestCMBGibbs;
    else if(name == "like_high")
        test = new TestLikeHigh;
    else if(name == "like_low")
        test = new TestLikeLow;
#endif
#endif
#endif
#ifdef COSMO_MINUIT
    else if(name == "fit")
        test = new TestFit;
#endif
#ifdef COSMO_CLASS
#ifdef COSMO_PLANCK
    else if(name == "planck_like")
        test = new TestPlanckLike;
#endif
#endif
#ifdef COSMO_CLASS
#ifdef COSMO_WMAP9
    else if(name == "wmap9_like")
        test = new TestWMAP9Like;
#endif
#endif
    else if(name == "wigner_3j")
        test = new TestWigner3J;
    else if(name == "table_function")
        test = new TestTableFunction;
    else if(name == "cubic_spline")
        test = new TestCubicSpline;
    else if(name == "three_rotation")
        test = new TestThreeRotation;
#ifdef COSMO_HEALPIX
    else if(name == "mask_apodizer")
        test = new TestMaskApodizer(1e-3);
#endif
#ifdef COSMO_CLASS
#ifdef COSMO_LAPACKPP
    else if(name == "matter_likelihood")
        test = new TestMatterLikelihood(1e-3);
#endif
#endif
#ifdef COSMO_ANN
    else if(name == "k_nearest_neighbors")
        test = new TestKNearestNeighbors;
#endif

    return test;
}

bool runTest(const std::string& name)
{
    TestFramework* test = createTest(name);
    check(test, "The test name was not found");
    unsigned int pass, fail;
    bool res = test->run(pass, fail);
    delete test;
    return res;
}

int main(int argc, char *argv[])
{
    try {
        StandardException exc;

        bool isMaster = CosmoMPI::create().isMaster();
        
        unsigned int total = 0, pass = 0, fail = 0;

        if(isMaster && argc < 2)
        {
            std::string exceptionStr = "The name of the test must be specified. Other options are 'all' to run all of the tests, 'fast' to run only the fast tests, 'slow' to run all of the slow tests, and 'list' to get a list of all the tests.";
            exc.set(exceptionStr);
            throw exc;
        }

        std::string argument(argv[1]);
        std::set<std::string> fastTests, slowTests;

        fastTests.insert("unit_conversions");
        fastTests.insert("int_operations");
        fastTests.insert("integral");
        fastTests.insert("conjugate_gradient");
        fastTests.insert("polynomial");
        fastTests.insert("legendre");
        fastTests.insert("spherical_harmonics");
#ifdef COSMO_LAPACKPP
        fastTests.insert("mcmc_fast");
#endif
#ifdef COSMO_MULTINEST
        fastTests.insert("multinest_fast");
#endif
#ifdef COSMO_CLASS
        fastTests.insert("cmb");
#endif
#ifdef COSMO_MINUIT
        fastTests.insert("fit");
#endif
#ifdef COSMO_CLASS
#ifdef COSMO_PLANCK
        fastTests.insert("planck_like");
#endif
#endif
#ifdef COSMO_CLASS
#ifdef COSMO_WMAP9
        fastTests.insert("wmap9_like");
#endif
#endif
        fastTests.insert("wigner_3j");
        fastTests.insert("table_function");
        fastTests.insert("cubic_spline");
        fastTests.insert("three_rotation");
#ifdef COSMO_CLASS
#ifdef COSMO_LAPACKPP
        fastTests.insert("matter_likelihood");
#endif
#endif
#ifdef COSMO_ANN
        fastTests.insert("k_nearest_neighbors");
#endif

#ifdef COSMO_PLANCK
#ifdef COSMO_CLASS
#ifdef COSMO_LAPACKPP
        slowTests.insert("mcmc_planck");
#endif
#ifdef COSMO_MULTINEST
        slowTests.insert("multinest_planck");
#endif
#endif
#endif
#ifdef COSMO_HEALPIX
#ifdef COSMO_LAPACKPP
#ifdef COSMO_CLASS
        slowTests.insert("cmb_gibbs");
        slowTests.insert("like_high");
        slowTests.insert("like_low");
#endif
#endif
#endif

#ifdef COSMO_HEALPIX
        slowTests.insert("mask_apodizer");
#endif

        if(argument == "all")
        {
            for(std::set<std::string>::const_iterator it = fastTests.begin(); it != fastTests.end(); ++it)
            {
                ++total;
                if(runTest(*it))
                    ++pass;
                else
                    ++fail;
            }
            for(std::set<std::string>::const_iterator it = slowTests.begin(); it != slowTests.end(); ++it)
            {
                ++total;
                if(runTest(*it))
                    ++pass;
                else
                    ++fail;
            }
        }
        else if(argument == "fast")
        {
            for(std::set<std::string>::const_iterator it = fastTests.begin(); it != fastTests.end(); ++it)
            {
                ++total;
                if(runTest(*it))
                    ++pass;
                else
                    ++fail;
            }
        }
        else if(argument == "slow")
        {
            for(std::set<std::string>::const_iterator it = slowTests.begin(); it != slowTests.end(); ++it)
            {
                ++total;
                if(runTest(*it))
                    ++pass;
                else
                    ++fail;
            }
        }
        else if(argument == "list")
        {
            if(isMaster)
            {
                output_screen_clean(std::endl << "FAST TESTS:" << std::endl);
                for(std::set<std::string>::const_iterator it = fastTests.begin(); it != fastTests.end(); ++it)
                    output_screen_clean("   " << *it << std::endl);

                output_screen_clean(std::endl << "SLOW TESTS:" << std::endl);
                for(std::set<std::string>::const_iterator it = slowTests.begin(); it != slowTests.end(); ++it)
                    output_screen_clean("   " << *it << std::endl);
                output_screen_clean(std::endl);
            }

            return 0;
        }
        else
        {
            bool found = false;
            std::set<std::string>::const_iterator it = fastTests.find(argument);
            if(it != fastTests.end())
                found = true;
            if(!found)
            {
                std::set<std::string>::const_iterator it = slowTests.find(argument);
                if(it != slowTests.end())
                    found = true;
            }

            if(!found)
            {
                if(isMaster)
                {
                    std::cout << "The test name " << argument << " was not found!" << std::endl;
                    std::cout << "Try 'all' to run all of the tests, 'fast' to run only the fast tests, 'slow' to run all of the slow tests, and 'list' to get a list of all the tests." << std::endl;
                }
                return -1;
            }

            ++total;
            if(runTest(argument))
                ++pass;
            else
                ++fail;
        }

        if(isMaster)
        {
            output_screen_clean(std::endl << "TOTAL NUMBER OF TESTS RUN: " << total << std::endl);
            output_screen_clean("PASSES: " << pass << std::endl);
            output_screen_clean("FAILURES: " << fail << std::endl);

            if(fail == 0)
            {
                output_screen_clean(std::endl << "\033[1;32mSUCCESS\033[0m" << std::endl << std::endl);
            }
            else
            {
                output_screen_clean(std::endl << "\033[1;31mFAIL\033[0m" << std::endl << std::endl);
            }
        }

    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}

