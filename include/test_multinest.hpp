#ifndef COSMO_PP_TEST_MULTINEST_HPP
#define COSMO_PP_TEST_MULTINEST_HPP

#include <test_framework.hpp>

class TestMultinestFast : public TestFramework
{
public:
    ~TestMultinestFast() {}

protected:
    bool isParallel(unsigned int i) const { return true; }
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

