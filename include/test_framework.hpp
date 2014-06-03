#ifndef COSMO_PP_TEST_FRAMEWORK_HPP
#define COSMO_PP_TEST_FRAMEWORK_HPP

#include <string>

class TestFramework
{
public:
    TestFramework(double precision = 1e-5);

    ~TestFramework() {}

    bool run(unsigned int& pass, unsigned int& fail);

    bool isMaster() const;

protected:
    virtual bool isParallel(unsigned int i) const { return false; }
    virtual std::string name() const = 0;
    virtual unsigned int numberOfSubtests() const = 0;
    virtual void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName) = 0;

protected:
    const double precision_;
};

#endif

