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
    void runSubTest10(double& res, double& expected, std::string& subTestName);
    void runSubTest11(double& res, double& expected, std::string& subTestName);
    void runSubTest12(double& res, double& expected, std::string& subTestName);
    void runSubTest13(double& res, double& expected, std::string& subTestName);
    void runSubTest14(double& res, double& expected, std::string& subTestName);
    void runSubTest15(double& res, double& expected, std::string& subTestName);
    void runSubTest16(double& res, double& expected, std::string& subTestName);
    void runSubTest17(double& res, double& expected, std::string& subTestName);
    void runSubTest18(double& res, double& expected, std::string& subTestName);
    void runSubTest19(double& res, double& expected, std::string& subTestName);

    void runSubTestEigen(double& res, double& expected, std::string& subTestName, bool pd);
};

#endif

