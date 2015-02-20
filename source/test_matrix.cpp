#include <macros.hpp>
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
    return 22;
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
    case 10:
        runSubTest10(res, expected, subTestName);
        break;
    case 11:
        runSubTest11(res, expected, subTestName);
        break;
    case 12:
        runSubTest12(res, expected, subTestName);
        break;
    case 13:
        runSubTest13(res, expected, subTestName);
        break;
    case 14:
        runSubTest14(res, expected, subTestName);
        break;
    case 15:
        runSubTest15(res, expected, subTestName);
        break;
    case 16:
        runSubTest16(res, expected, subTestName);
        break;
    case 17:
        runSubTest17(res, expected, subTestName);
        break;
    case 18:
        runSubTest18(res, expected, subTestName);
        break;
    case 19:
        runSubTest19(res, expected, subTestName);
        break;
    case 20:
        runSubTestEigen(res, expected, subTestName, false);
        break;
    case 21:
        runSubTestEigen(res, expected, subTestName, true);
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
    subTestName = "simple_subtract";
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

void
TestMatrix::runSubTest10(double& res, double& expected, std::string& subTestName)
{
    Math::SymmetricMatrix<int> mat(3, -1);

    mat(2, 1) = 7;
    res = mat(1, 2);
    expected = mat(2, 1);
    subTestName = "simple_symmetric";
}

void
TestMatrix::runSubTest11(double& res, double& expected, std::string& subTestName)
{
    Math::SymmetricMatrix<int> mat(3);

    mat(0, 1) = -7;

    std::string fileName = "test_files/matrix_test_11.dat";

    mat.writeIntoFile(fileName.c_str());

    Math::SymmetricMatrix<int> mat2(fileName.c_str());

    res = mat2(1, 0);
    expected = mat(0, 1);
    subTestName = "simple_symmetric_read_write";
}

void
TestMatrix::runSubTest12(double& res, double& expected, std::string& subTestName)
{
    Math::SymmetricMatrix<int> mat(3);

    mat(0, 1) = -7;

    std::string fileName = "test_files/matrix_test_12.txt";

    mat.writeIntoTextFile(fileName.c_str());

    Math::SymmetricMatrix<int> mat2;
    mat2.readFromTextFile(fileName.c_str());

    res = mat2(1, 0);
    expected = mat(0, 1);
    subTestName = "simple_symmetric_read_write_text";
}

void
TestMatrix::runSubTest13(double& res, double& expected, std::string& subTestName)
{
    Math::SymmetricMatrix<int> mat(4);

    mat(0, 1) = -7;

    Math::SymmetricMatrix<int> mat2;
    mat2 = mat;

    res = mat2(1, 0);
    expected = mat(0, 1);
    subTestName = "simple_symmetric_copy";
}

void
TestMatrix::runSubTest14(double& res, double& expected, std::string& subTestName)
{
    Math::SymmetricMatrix<int> mat(2, 5), mat1(2, 3);

    mat(0, 1) = -7;
    mat1(1, 0) = 5;

    Math::SymmetricMatrix<int> mat2 = mat + mat1;

    res = mat2(1, 0);
    expected = mat(1, 0) + mat1(0, 1);
    subTestName = "simple_symmetric_add";
}

void
TestMatrix::runSubTest15(double& res, double& expected, std::string& subTestName)
{
    Math::SymmetricMatrix<int> mat(2, 10), mat1(2, -20);

    mat(0, 1) = -7;
    mat1(0, 1) = 5;

    Math::SymmetricMatrix<int> mat2 = mat;
    mat2 -= mat1;

    res = mat2(0, 1);
    expected = mat(1, 0) - mat1(0, 1);
    subTestName = "simple_symmetric_subtract";
}

void
TestMatrix::runSubTest16(double& res, double& expected, std::string& subTestName)
{
    Math::SymmetricMatrix<int> mat(3);

    mat(1, 0) = -7;

    Math::Matrix<int> mat2 = mat;

    res = mat2(0, 1);
    expected = mat(0, 1);

    subTestName = "simple_matrix_from_symmetric";
}

void
TestMatrix::runSubTest17(double& res, double& expected, std::string& subTestName)
{
    Math::SymmetricMatrix<double> mat(3, 5), mat1(3, 7);

    Math::SymmetricMatrix<double> mat2 = mat * mat1;

    res = mat2(1, 2);
    expected = mat(0, 1) * mat1(2, 0) * mat.size();
    subTestName = "simple_symmetric_multiply";
}

void
TestMatrix::runSubTest18(double& res, double& expected, std::string& subTestName)
{
    expected = 1;
    subTestName = "simple_symmetric_invert";

#ifdef COSMO_LAPACK
    Math::SymmetricMatrix<double> mat(2);
    mat(0, 0) = 2;
    mat(1, 1) = 3;
    mat(1, 0) = 1;

    mat.writeIntoTextFile("test_files/matrix_test_18_original.txt");

    Math::SymmetricMatrix<double> invMat = mat;
    invMat.invert();

    invMat.writeIntoTextFile("test_files/matrix_test_18_inverse.txt");

    Math::SymmetricMatrix<double> prod = mat;
    prod *= invMat;
    prod.writeIntoTextFile("test_files/matrix_test_18_product.txt");

    res = 1;
    for(int i = 0; i < prod.size(); ++i)
    {
        for(int j = 0; j < prod.size(); ++j)
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
TestMatrix::runSubTest19(double& res, double& expected, std::string& subTestName)
{
    subTestName = "simple_symmetric_determinant";

#ifdef COSMO_LAPACK
    Math::SymmetricMatrix<double> mat(2);
    mat(0, 0) = 3;
    mat(1, 1) = 4;
    mat(0, 1) = 2;

    res = mat.determinant();
    expected = mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0);

#else
    output_screen_clean("This test (below) is skipped because Cosmo++ has not been linked to lapack" << std::endl);
    res = 1;
    expected = 1;
#endif
}

void
TestMatrix::runSubTestEigen(double& res, double& expected, std::string& subTestName, bool pd)
{
    subTestName = "simple_symmetric_eigen";
    if(pd)
        subTestName = "simple_symmetric_positive_eigen";

#ifdef COSMO_LAPACK
    Math::SymmetricMatrix<double> mat(3);
    mat(0, 0) = 3;
    mat(1, 1) = 10;
    mat(2, 2) = 4;
    mat(0, 2) = 2;

    std::vector<double> eigenvals;
    Math::Matrix<double> eigenvecs;

    const int info = mat.getEigen(&eigenvals, &eigenvecs, pd);
    if(info)
    {
        output_screen_clean("FAIL! Eigenvalue/eigenvector decomposition failed. Info = " << info << std::endl);
        res = 0;
        return;
    }

    //output_screen_clean("Eigenvalues: " << eigenvals[0] << ", " << eigenvals[1] << ", " << eigenvals[2] << std::endl);

    if(pd)
        eigenvecs.writeIntoTextFile("test_files/matrix_test_eigenvecs.txt");
    else
        eigenvecs.writeIntoTextFile("test_files/matrix_test_eigenvecs_pos.txt");

    res = 1;
    expected = 1;

    Math::Matrix<double> m = mat;

    for(int i = 0; i < 3; ++i)
    {
        Math::Matrix<double> v = eigenvecs.getCol(i);
        Math::Matrix<double> prod = m * v;
        for(int j = 0; j < 3; ++j)
        {
            if(!Math::areEqual(prod(j, 0), eigenvals[i] * v(j, 0), 1e-5))
            {
                output_screen_clean("FAIL! The eigenvalue " << i << " times the eigenvector doesn't match the matrix times the eigenvector." << std::endl);
                output_screen_clean("\tLooking at index " << j << ", expected " << eigenvals[i] * v(j, 0) << " obtained " << prod(j, 0) << std::endl);
                res = 0;
            }
        }
    }

    Math::Matrix<double> diag = eigenvecs.getTranspose() * mat * eigenvecs;

    for(int i = 0; i < 3; ++i)
    {
        for(int j = 0; j < 3; ++j)
        {
            if(i == j)
            {
                if(!Math::areEqual(diag(i, i), eigenvals[i], 1e-5))
                {
                    output_screen_clean("FAIL! The diagonalized matrix has " << diag(i, i) << " on the diagonal at index " << i << " but the corresponding eigenvalue is " << eigenvals[i] << std::endl);
                    res = 0;
                }
            }
            else
            {
                if(!Math::areEqual(diag(i, j), 0.0, 1e-5))
                {
                    output_screen_clean("FAIL! The diagonalized matrix has " << diag(i, j) << " as the off-diagonal element (" << i << ", " << j << "). Must be 0." << std::endl);
                    res = 0;
                }
            }
        }
    }

#else
    output_screen_clean("This test (below) is skipped because Cosmo++ has not been linked to lapack" << std::endl);
    res = 1;
    expected = 1;
#endif
}

