#ifndef COSMO_PP_TEST_UNIT_CONVERSIONS_HPP
#define COSMO_PP_TEST_UNIT_CONVERSIONS_HPP

#include <test_framework.hpp>

class TestUnitConversions : public TestFramework
{
public:
    ~TestUnitConversions() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

