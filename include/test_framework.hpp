#ifndef COSMO_PP_TEST_FRAMEWORK_HPP
#define COSMO_PP_TEST_FRAMEWORK_HPP

#include <string>

class TestFramework
{
public:
    /// Constructor.
    /// \param precision The precision with which the actual and expected numbers should be compared.
    TestFramework(double precision = 1e-5);

    /// Destructor.
    ~TestFramework() {}

    /// Run all of the subtests.
    /// \param pass The number of passed subtests upon return.
    /// \param fail The number of failed subtests upon return.
    /// \return true if all of the subtests passed.
    bool run(unsigned int& pass, unsigned int& fail);

    /// Check if this is the master MPI process.
    bool isMaster() const;

protected:
    /// Tells if the test should be run in parallel or not. If false, when running on MPI only the master process will run it.
    virtual bool isParallel(unsigned int i) const { return false; }

    /// The name of the test. Purely virtual.
    virtual std::string name() const = 0;

    /// The number of subtests. Purely virtual.
    virtual unsigned int numberOfSubtests() const = 0;

    /// Run a given subtest. Purely virtual.
    /// \param i The index of the subtest, starting from 0.
    /// \param res The result of the subtest to be compared to the expected result.
    /// \param expected The expected result of the subtest.
    /// \param subTestName The name of the subtest.
    virtual void runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName) = 0;

protected:
    const double precision_;
};

#endif

