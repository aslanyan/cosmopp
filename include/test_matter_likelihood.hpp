#ifndef COSMO_PP_TEST_MATTER_LIKELIHOOD_HPP
#define COSMO_PP_TEST_MATTER_LIKELIHOOD_HPP

#include <test_framework.hpp>

class TestMatterLikelihood : public TestFramework
{
public:
    ~TestMatterLikelihood() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

