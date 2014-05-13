#ifndef COSMO_PP_TEST_INTEGRAL_HPP
#define COSMO_PP_TEST_INTEGRAL_HPP

#include <test_framework.hpp>

class TestIntegral : public TestFramework
{
public:
    ~TestIntegral() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

