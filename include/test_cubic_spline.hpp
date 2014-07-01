#ifndef COSMO_PP_TEST_CUBIC_SPLINE_HPP
#define COSMO_PP_TEST_CUBIC_SPLINE_HPP

#include <test_framework.hpp>

class TestCubicSpline : public TestFramework
{
public:
    ~TestCubicSpline() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif


