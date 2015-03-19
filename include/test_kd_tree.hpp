#ifndef COSMO_PP_TEST_KD_TREE_HPP
#define COSMO_PP_TEST_KD_TREE_HPP

#include <test_framework.hpp>

class TestKDTree : public TestFramework
{
public:
    ~TestKDTree() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
    void runSubTest0(double& res, double& expected, std::string& subTestName);
    void runSubTest1(double& res, double& expected, std::string& subTestName);
    void runSubTest2(double& res, double& expected, std::string& subTestName);
    void runSubTest3(double& res, double& expected, std::string& subTestName);
    void runSubTest4(double& res, double& expected, std::string& subTestName);
};

#endif

