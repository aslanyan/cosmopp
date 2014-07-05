#ifndef COSMO_PP_TEST_LIKE_LOW_HPP
#define COSMO_PP_TEST_LIKE_LOW_HPP

#include <test_framework.hpp>

class TestLikeLow : public TestFramework
{
public:
    ~TestLikeLow() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

