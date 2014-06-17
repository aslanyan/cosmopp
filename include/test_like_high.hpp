#ifndef COSMO_PP_TEST_LIKE_HIGH_HPP
#define COSMO_PP_TEST_LIKE_HIGH_HPP

#include <test_framework.hpp>

class TestLikeHigh : public TestFramework
{
public:
    ~TestLikeHigh() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

