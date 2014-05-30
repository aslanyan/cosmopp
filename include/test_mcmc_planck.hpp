#ifndef COSMO_PP_TEST_MCMC_PLANCK_HPP
#define COSMO_PP_TEST_MCMC_PLANCK_HPP

#include <test_framework.hpp>

class TestMCMCPlanck : public TestFramework
{
public:
    ~TestMCMCPlanck() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

