#ifndef COSMO_PP_TEST_CMB_HPP
#define COSMO_PP_TEST_CMB_HPP

#include <test_framework.hpp>

class TestCMB : public TestFramework
{
public:
    ~TestCMB() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

