#ifndef COSMO_PP_TEST_LEGENDRE_HPP
#define COSMO_PP_TEST_LEGENDRE_HPP

#include <test_framework.hpp>

class TestLegendre : public TestFramework
{
public:
    ~TestLegendre() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

