#ifndef COSMO_PP_TEST_MULTINEST_PLANCK_HPP
#define COSMO_PP_TEST_MULTINEST_PLANCK_HPP

#include <test_framework.hpp>

class TestMultinestPlanck : public TestFramework
{
public:
    ~TestMultinestPlanck() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

