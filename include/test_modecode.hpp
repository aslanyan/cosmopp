#ifndef COSMO_PP_TEST_MODECODE_HPP
#define COSMO_PP_TEST_MODECODE_HPP

#include <test_framework.hpp>

class TestModeCode : public TestFramework
{
public:
    ~TestModeCode() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

