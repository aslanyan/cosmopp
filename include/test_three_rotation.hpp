#ifndef COSMO_PP_TEST_THREE_ROTATION_HPP
#define COSMO_PP_TEST_THREE_ROTATION_HPP

#include <test_framework.hpp>

class TestThreeRotation : public TestFramework
{
public:
    ~TestThreeRotation() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

