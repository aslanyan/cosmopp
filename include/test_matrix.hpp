#ifndef COSMO_PP_TEST_MATRIX_HPP
#define COSMO_PP_TEST_MATRIX_HPP

#include <test_framework.hpp>

class TestMatrix : public TestFramework
{
public:
    ~TestMatrix() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

