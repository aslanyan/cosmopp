#ifndef COSMO_PP_TEST_FAST_APPROXIMATOR_ERROR_HPP
#define COSMO_PP_TEST_FAST_APPROXIMATOR_ERROR_HPP

#include <test_framework.hpp>

class TestFastApproximatorError : public TestFramework
{
public:
    TestFastApproximatorError(double precision = 1e-3) : TestFramework(precision) {}
    ~TestFastApproximatorError() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

