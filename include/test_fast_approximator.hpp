#ifndef COSMO_PP_TEST_FAST_APPROXIMATOR_HPP
#define COSMO_PP_TEST_FAST_APPROXIMATOR_HPP

#include <test_framework.hpp>

class TestFastApproximator : public TestFramework
{
public:
    TestFastApproximator(double precision = 1e-3) : TestFramework(precision) {}
    ~TestFastApproximator() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

