#ifndef COSMO_PP_TEST_SPHERICAL_HARMONICS_HPP
#define COSMO_PP_TEST_SPHERICAL_HARMONICS_HPP

#include <test_framework.hpp>

class TestSphericalHarmonics : public TestFramework
{
public:
    ~TestSphericalHarmonics() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

