#ifndef COSMO_PP_TEST_WMAP9_LIKE_HPP
#define COSMO_PP_TEST_WMAP9_LIKE_HPP

#include <test_framework.hpp>

class TestWMAP9Like : public TestFramework
{
public:
    ~TestWMAP9Like() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif
