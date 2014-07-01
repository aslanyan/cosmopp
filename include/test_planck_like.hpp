#ifndef COSMO_PP_TEST_PLANCK_LIKE_HPP
#define COSMO_PP_TEST_PLANCK_LIKE_HPP

#include <test_framework.hpp>

class TestPlanckLike : public TestFramework
{
public:
    ~TestPlanckLike() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

