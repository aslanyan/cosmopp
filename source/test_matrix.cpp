#include <matrix.hpp>
#include <test_matrix.hpp>

std::string
TestMatrix::name() const
{
    return std::string("MATRIX TESTER");
}

unsigned int
TestMatrix::numberOfSubtests() const
{
    return 1;
}

void
TestMatrix::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    Math::Matrix<int> mat(2, 3, -1);

    res = mat(1, 2);
    expected = -1;
    subTestName = "simple1";
}
