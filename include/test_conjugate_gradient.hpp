#ifndef COSMO_PP_TEST_CONJUGATE_GRADIENT_HPP
#define COSMO_PP_TEST_CONJUGATE_GRADIENT_HPP
#include <test_framework.hpp>

class TestConjugateGradient : public TestFramework
{
public:
    ~TestConjugateGradient() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

