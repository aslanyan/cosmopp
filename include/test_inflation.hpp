#ifndef COSMO_PP_TEST_INFLATION_HPP
#define COSMO_PP_TEST_INFLATION_HPP

#include <test_framework.hpp>

class TestInflation : public TestFramework
{
public:
    ~TestInflation() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

