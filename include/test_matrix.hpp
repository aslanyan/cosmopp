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

private:
    void runSubTest0(double& res, double& expected, std::string& subTestName);
    void runSubTest1(double& res, double& expected, std::string& subTestName);
    void runSubTest2(double& res, double& expected, std::string& subTestName);
    void runSubTest3(double& res, double& expected, std::string& subTestName);
    void runSubTest4(double& res, double& expected, std::string& subTestName);
    void runSubTest5(double& res, double& expected, std::string& subTestName);
    void runSubTest6(double& res, double& expected, std::string& subTestName);
    void runSubTest7(double& res, double& expected, std::string& subTestName);
    void runSubTest8(double& res, double& expected, std::string& subTestName);
    void runSubTest9(double& res, double& expected, std::string& subTestName);
};

#endif

