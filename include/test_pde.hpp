#ifndef COSMO_PP_TEST_PDE_HPP
#define COSMO_PP_TEST_PDE_HPP

#include <ctime>

#include <test_framework.hpp>
#include <random.hpp>

class TestPDE : public TestFramework
{
public:
    TestPDE(double precision = 1e-3) : TestFramework(precision), rand_(std::time(0), 0, 100000) {}
    ~TestPDE() {}

protected:
    std::string name() const;
    virtual bool isParallel(unsigned int i) const { return true; }
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);

private:
    Math::UniformIntGenerator rand_;
};

#endif

