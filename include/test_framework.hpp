#ifndef COSMO_PP_TEST_FRAMEWORK_HPP
#define COSMO_PP_TEST_FRAMEWORK_HPP

#include <string>

#include <macros.hpp>
#include <numerics.hpp>

class TestFramework
{
public:
    TestFramework(double precision = 1e-5) : precision_(precision)
    {
        check(precision > 0, "invalid precision " << precision << ", must be positive");
        check(precision < 1, "precision " << precision << " is too big, needs to be smaller than 1");
    }

    ~TestFramework()
    {
    }

    bool run(unsigned int& pass, unsigned int& fail)
    {
        output_screen(std::endl << "TEST: " << name() << std::endl << std::endl);
        pass = 0;
        fail = 0;
        const unsigned int n = numberOfSubtests();
        check(n > 0, "need at least 1 subtest");
        for(unsigned int i = 0; i < n; ++i)
        {
            double res, expected;
            std::string subTestName;
            runSubTest(i, res, expected, subTestName);

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

        output_screen(std::endl << "Total subtests: " << pass + fail << std::endl << "Pass: " << pass << std::endl << "Fail: " << fail << std::endl << std::endl);
        return (fail == 0);
    }

protected:
    virtual std::string name() const = 0;
    virtual unsigned int numberOfSubtests() const = 0;
    virtual void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName) = 0;

protected:
    const double precision_;
};

#endif

