#ifndef COSMO_PP_TEST_K_NEAREST_NEIGHBORS_HPP
#define COSMO_PP_TEST_K_NEAREST_NEIGHBORS_HPP

#include <test_framework.hpp>

class TestKNearestNeighbors : public TestFramework
{
public:
    ~TestKNearestNeighbors() {}

protected:
    std::string name() const;
    unsigned int numberOfSubtests() const;
    void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName);
};

#endif
