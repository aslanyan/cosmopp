#include <matrix_impl.hpp>
#include <test_matrix.hpp>
#include <numerics.hpp>

std::string
TestMatrix::name() const
{
    return std::string("MATRIX TESTER");
}

unsigned int
TestMatrix::numberOfSubtests() const
{
    return 10;
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
        break;
    case 2:
        runSubTest2(res, expected, subTestName);
        break;
    case 3:
        runSubTest3(res, expected, subTestName);
        break;
    case 4:
        runSubTest4(res, expected, subTestName);
        break;
    case 5:
        runSubTest5(res, expected, subTestName);
        break;
    case 6:
        runSubTest6(res, expected, subTestName);
        break;
    case 7:
        runSubTest7(res, expected, subTestName);
        break;
    case 8:
        runSubTest8(res, expected, subTestName);
        break;
    case 9:
        runSubTest9(res, expected, subTestName);
        break;
    default:
        check(false, "");
        break;
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

void
TestMatrix::runSubTest2(double& res, double& expected, std::string& subTestName)
{
    Math::Matrix<int> mat(2, 3);

    mat(0, 1) = -7;

    std::string fileName = "test_files/matrix_test_2.txt";

    mat.writeIntoTextFile(fileName.c_str());

    Math::Matrix<int> mat2;
    mat2.readFromTextFile(fileName.c_str());

    res = mat2(0, 1);
    expected = mat(0, 1);
    subTestName = "simple_read_write_text";
}

void
TestMatrix::runSubTest3(double& res, double& expected, std::string& subTestName)
{
    Math::Matrix<int> mat(2, 3);

    mat(0, 1) = -7;

    Math::Matrix<int> mat2;
    mat2 = mat;

    res = mat2(0, 1);
    expected = mat(0, 1);
    subTestName = "simple_copy";
}

void
TestMatrix::runSubTest4(double& res, double& expected, std::string& subTestName)
{
    Math::Matrix<int> mat(2, 3), mat1(2, 3);

    mat(0, 1) = -7;
    mat1(0, 1) = 5;

    Math::Matrix<int> mat2 = mat + mat1;

    res = mat2(0, 1);
    expected = mat(0, 1) + mat1(0, 1);
    subTestName = "simple_add";
}

void
TestMatrix::runSubTest5(double& res, double& expected, std::string& subTestName)
{
    Math::Matrix<int> mat(2, 3), mat1(2, 3);

    mat(0, 1) = -7;
    mat1(0, 1) = 5;

    Math::Matrix<int> mat2 = mat;
    mat2 -= mat1;

    res = mat2(0, 1);
    expected = mat(0, 1) - mat1(0, 1);
    subTestName = "simple_add";
}

void
TestMatrix::runSubTest6(double& res, double& expected, std::string& subTestName)
{
    Math::Matrix<int> mat(2, 3);

    mat(0, 1) = -7;

    Math::Matrix<int> mat2 = mat;
    mat2.transpose();

    res = mat2(1, 0);
    expected = mat(0, 1);
    subTestName = "simple_transpose";
}

void
TestMatrix::runSubTest7(double& res, double& expected, std::string& subTestName)
{
    Math::Matrix<double> mat(2, 3, 5), mat1(3, 4, 7);

    Math::Matrix<double> mat2 = mat * mat1;

    res = mat2(1, 3);
    expected = mat(0, 0) * mat1(0, 0) * mat.cols();
    subTestName = "simple_multiply";
}

void
TestMatrix::runSubTest8(double& res, double& expected, std::string& subTestName)
{
    expected = 1;
    subTestName = "simple_invert";

#ifdef COSMO_LAPACK
    Math::Matrix<double> mat(2, 2);
    mat(0, 0) = 1;
    mat(0, 1) = 2;
    mat(1, 0) = 3;
    mat(1, 1) = 4;

    mat.writeIntoTextFile("test_files/matrix_test_8_original.txt");

    Math::Matrix<double> invMat = mat.getInverse();

    invMat.writeIntoTextFile("test_files/matrix_test_8_inverse.txt");

    Math::Matrix<double> prod = invMat * mat;
    prod.writeIntoTextFile("test_files/matrix_test_8_product.txt");

    res = 1;
    for(int i = 0; i < prod.rows(); ++i)
    {
        for(int j = 0; j < prod.cols(); ++j)
        {
            if(i == j)
            {
                if(!Math::areEqual(prod(i, j), 1.0, 1e-5))
                {
                    output_screen("FAIL! Diagonal element " << i << " must be 1 but it is " << prod(i, j) << std::endl);
                    res = 0;
                }
            }
            else
            {
                if(!Math::areEqual(prod(i, j), 0.0, 1e-5))
                {
                    output_screen("FAIL! Non-diagonal element " << i << " " << j << " must be 0 but it is " << prod(i, j) << std::endl);
                    res = 0;
                }
            }
        }
    }
#else
    output_screen_clean("This test (below) is skipped because Cosmo++ has not been linked to lapack" << std::endl);
    res = 1;
#endif
}

void
TestMatrix::runSubTest9(double& res, double& expected, std::string& subTestName)
{
    subTestName = "simple_determinant";

#ifdef COSMO_LAPACK
    Math::Matrix<double> mat(2, 2);
    mat(0, 0) = 1;
    mat(0, 1) = 2;
    mat(1, 0) = 3;
    mat(1, 1) = 4;

    res = mat.determinant();
    expected = mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0);

#else
    output_screen_clean("This test (below) is skipped because Cosmo++ has not been linked to lapack" << std::endl);
    res = 1;
    expected = 1;
#endif
}

