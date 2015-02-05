#ifndef COSMO_PP_TEST_GAUSSIAN_PROCESS_HPP
#define COSMO_PP_TEST_GAUSSIAN_PROCESS_HPP

#include <test_framework.hpp>

class TestGaussianProcess : public TestFramework
{
public:
    ~TestGaussianProcess() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);

private:
    void runSubTest0(double& res, double& expected, std::string& subTestName);
    void runSubTest1(double& res, double& expected, std::string& subTestName);
};

#endif

