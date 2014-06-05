#ifndef COSMO_PP_TEST_CMB_GIBBS_HPP
#define COSMO_PP_TEST_CMB_GIBBS_HPP

#include <test_framework.hpp>

class TestCMBGibbs : public TestFramework
{
public:
    ~TestCMBGibbs() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

