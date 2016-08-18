#include <cosmo_mpi.hpp>
#include <macros.hpp>
#include <numerics.hpp>
#include <test_framework.hpp>

TestFramework::TestFramework(double precision) : precision_(precision)
{
    check(precision > 0, "invalid precision " << precision << ", must be positive");
    check(precision < 1, "precision " << precision << " is too big, needs to be smaller than 1");
}

bool
TestFramework::run(unsigned int& pass, unsigned int& fail)
{
    CosmoMPI::create().barrier();
    if(isMaster())
    {
        output_screen_clean(std::endl << "TEST: " << name() << std::endl << std::endl);
        pass = 0;
        fail = 0;
    }
    const unsigned int n = numberOfSubtests();
    check(n > 0, "need at least 1 subtest");
    for(unsigned int i = 0; i < n; ++i)
    {
        double res, expected;
        std::string subTestName;
        if(isMaster() || isParallel(i))
            runSubTest(i, res, expected, subTestName);

        if(isMaster())
        {
            const bool testRes = Math::areEqual(expected, res, precision_);
            output_screen_clean("   " << subTestName << ": ");
            if(testRes)
            {
                ++pass;
                output_screen_clean("\033[1;32mPASS\033[0m" << std::endl);
            }
            else
            {
                ++fail;
                output_screen_clean("\033[1;31mFAIL\033[0m" << std::endl);
                output_screen_clean("      Result: " << res << "   Expected result: " << expected << std::endl);
            }
        }
    }

    CosmoMPI::create().barrier();
    if(isMaster())
    {
        output_screen_clean(std::endl << "Total subtests: " << pass + fail << std::endl << "Pass: " << pass << std::endl << "Fail: " << fail << std::endl << std::endl);
        return (fail == 0);
    }

    return true;
}

bool
TestFramework::isMaster() const
{
    return CosmoMPI::create().isMaster();
}
