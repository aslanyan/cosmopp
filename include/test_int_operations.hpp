#ifndef COSMO_PP_TEST_INT_OPERATIONS_HPP
#define COSMO_PP_TEST_INT_OPERATIONS_HPP

#include <test_framework.hpp>

class TestIntOperations : public TestFramework
{
public:
    ~TestIntOperations() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

