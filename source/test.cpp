#include <fstream>
#include <set>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>

#include <test_framework.hpp>
#include <test_unit_conversions.hpp>
#include <test_int_operations.hpp>
#include <test_integral.hpp>
#include <test_conjugate_gradient.hpp>
#include <test_polynomial.hpp>
#include <test_legendre.hpp>
#include <test_mcmc.hpp>

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
    else if(name == "mcmc_fast")
        test = new TestMCMCFast;

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
        
        unsigned int total = 0, pass = 0, fail = 0;

        if(argc < 2)
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
        fastTests.insert("mcmc_fast");

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
            std::cout << std::endl << "FAST TESTS:" << std::endl;
            for(std::set<std::string>::const_iterator it = fastTests.begin(); it != fastTests.end(); ++it)
                std::cout << "   " << *it << std::endl;

            std::cout << std::endl << "SLOW TESTS:" << std::endl;
            for(std::set<std::string>::const_iterator it = slowTests.begin(); it != slowTests.end(); ++it)
                std::cout << "   " << *it << std::endl;
            std::cout << std::endl;
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
                std::cout << "The test name " << argument << " was not found!" << std::endl;
                std::cout << "Try 'all' to run all of the tests, 'fast' to run only the fast tests, 'slow' to run all of the slow tests, and 'list' to get a list of all the tests." << std::endl;
                return -1;
            }

            ++total;
            if(runTest(argument))
                ++pass;
            else
                ++fail;
        }

        std::cout << std::endl << "TOTAL NUMBER OF TESTS RUN: " << total << std::endl;
        std::cout << "PASSES: " << pass << std::endl;
        std::cout << "FAILURES: " << fail << std::endl;

        if(fail == 0)
            std::cout << std::endl << "\033[1;32mSUCCESS\033[0m" << std::endl << std::endl;
        else
            std::cout << std::endl << "\033[1;31mFAIL\033[0m" << std::endl << std::endl;
        
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}

