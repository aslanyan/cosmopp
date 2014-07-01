#ifndef COSMO_PP_TEST_TABLE_FUNCTION_HPP
#define COSMO_PP_TEST_TABLE_FUNCTION_HPP

#include <test_framework.hpp>

class TestTableFunction : public TestFramework
{
public:
    ~TestTableFunction() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif

