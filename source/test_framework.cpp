#ifdef COSMO_MPI
#include <mpi.h>
#endif

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
    if(isMaster())
    {
        output_screen(std::endl << "TEST: " << name() << std::endl << std::endl);
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
            output_screen("   " << subTestName << ": ");
            if(testRes)
            {
                ++pass;
                output_screen("\033[1;32mPASS\033[0m" << std::endl);
            }
            else
            {
                ++fail;
                output_screen("\033[1;31mFAIL\033[0m" << std::endl);
                output_screen("      Result: " << res << "   Expected result: " << expected << std::endl);
            }
        }
    }

    if(isMaster())
    {
        output_screen(std::endl << "Total subtests: " << pass + fail << std::endl << "Pass: " << pass << std::endl << "Fail: " << fail << std::endl << std::endl);
        return (fail == 0);
    }
    return true;
}

bool
TestFramework::isMaster() const
{
#ifdef COSMO_MPI
    int hasMpiInitialized;
    MPI_Initialized(&hasMpiInitialized);
    if(!hasMpiInitialized)
        MPI_Init(NULL, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return (rank == 0);
#endif

    return true;
}
