#ifndef COSMO_PP_TEST_MULTINEST_HPP
#define COSMO_PP_TEST_MULTINEST_HPP

#include <test_framework.hpp>

class TestMultinestFast : public TestFramework
{
public:
    ~TestMultinestFast() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

