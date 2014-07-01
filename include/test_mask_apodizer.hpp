#ifndef COSMO_PP_TEST_MASK_APODIZER_HPP
#define COSMO_PP_TEST_MASK_APODIZER_HPP

#include <test_framework.hpp>

class TestMaskApodizer : public TestFramework
{
public:
    TestMaskApodizer(double precision = 1e-5) : TestFramework(precision) {}

    ~TestMaskApodizer() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

