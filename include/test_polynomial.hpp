#ifndef COSMO_PP_TEST_POLYNOMIAL_HPP
#define COSMO_PP_TEST_POLYNOMIAL_HPP

#include <test_framework.hpp>

class TestPolynomial : public TestFramework
{
public:
    ~TestPolynomial() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

