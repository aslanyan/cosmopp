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
    return 2;
}

void
TestMatrix::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    switch(i)
    {
    case 0:
        runSubTest0(res, expected, subTestName);
        break;
    case 1:
        runSubTest1(res, expected, subTestName);
    }
}

void
TestMatrix::runSubTest0(double& res, double& expected, std::string& subTestName)
{
    Math::Matrix<int> mat(2, 3, -1);

    res = mat(1, 2);
    expected = -1;
    subTestName = "simple_constructor";
}

void
TestMatrix::runSubTest1(double& res, double& expected, std::string& subTestName)
{
    Math::Matrix<int> mat(2, 3);

    mat(0, 1) = -7;

    std::string fileName = "test_files/matrix_test_1.dat";

    mat.writeIntoFile(fileName.c_str());

    Math::Matrix<int> mat2(fileName.c_str());

    res = mat2(0, 1);
    expected = mat(0, 1);
    subTestName = "simple_read_write";
}
