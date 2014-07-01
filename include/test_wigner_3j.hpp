#ifndef COSMO_PP_TEST_WIGNER_3J_HPP
#define COSMO_PP_TEST_WIGNER_3J_HPP

#include <test_framework.hpp>

class TestWigner3J : public TestFramework
{
public:
    ~TestWigner3J() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

