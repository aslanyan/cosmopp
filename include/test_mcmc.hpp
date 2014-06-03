#ifndef COSMO_PP_TEST_MCMC_HPP
#define COSMO_PP_TEST_MCMC_HPP

#include <test_framework.hpp>

class TestMCMCFast : public TestFramework
{
public:
    ~TestMCMCFast() {}

protected:
    bool isParallel(unsigned int i) const { return true; }
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

