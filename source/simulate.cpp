#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <ctime>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <numerics.hpp>
#include <whole_matrix.hpp>
#include <simulate.hpp>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include <gmd.h>
#include <lavd.h>
#include <laslv.h>
#include <lavli.h>
#include <blas2pp.h>
#include <blas3pp.h>

#include "alm.h"
#include "xcomplex.h"
#include "chealpix.h"

int index(int l, int m, int lMin)
{
    return l * (l + 1) / 2 + m - lMin * (lMin + 1) / 2;
}

int index1(int l, int m, int lMin)
{
    return l * (l + 1) + m - (lMin * (lMin + 1) - lMin);
}

void
Simulate::simulateAlm(const WholeMatrix& wholeMatrix, Alm<xcomplex<double> >& alm, double* chi2, int* dof, time_t seed)
{
    if(seed == 0)
        seed = std::time(0);
    
    const int lMin = wholeMatrix.getLMin(), lMax = wholeMatrix.getLMax();
    
    
    
    //output_screen("Preparing the real and imaginary covariance matrices..." << std::endl);
    const int matrixSize = index(lMax, lMax, lMin) + 1;
    LaGenMatDouble reMatrix(matrixSize, matrixSize), imMatrix(matrixSize, matrixSize);
    
    for(int l1 = lMin; l1 <= lMax; ++l1)
    {
        for(int m1 = 0; m1 <= l1; ++m1)
        {
            for(int l = lMin; l <= lMax; ++l)
            {
                for(int m = 0; m <= l; ++m)
                {
                    const int i = index(l, m, lMin);
                    const int j = index(l1, m1, lMin);
                    check(i < matrixSize, "");
                    check(j < matrixSize, "");
                    const int minus1m = (m % 2 ? -1 : 1);
                    reMatrix(i, j) = (wholeMatrix.element(l1, m1, l, m) + minus1m * wholeMatrix.element(l1, m1, l, -m)) / 2;
                    imMatrix(i, j) = (wholeMatrix.element(l1, m1, l, m) - minus1m * wholeMatrix.element(l1, m1, l, -m)) / 2;
                }
            }
        }
    }
    //output_screen("OK" << std::endl);
    
    LaVectorDouble reEigenvalsRe(matrixSize), reEigenvalsIm(matrixSize), imEigenvalsRe(matrixSize), imEigenvalsIm(matrixSize);
    LaGenMatDouble reVecs(matrixSize, matrixSize), imVecs(matrixSize, matrixSize);

    diagonalizeMatrix(reMatrix, reEigenvalsRe, reEigenvalsIm, reVecs);
    diagonalizeMatrix(imMatrix, imEigenvalsRe, imEigenvalsIm, imVecs);
    
    //output_screen("Simulating with seed " << seed << "..." << std::endl)
    boost::variate_generator<boost::mt19937, boost::normal_distribution<> > generator(boost::mt19937(seed), boost::normal_distribution<>(0, 1));
    
    const xcomplex<double> zero(0, 0);
    LaVectorDouble re(matrixSize), im(matrixSize), reRot(matrixSize), imRot(matrixSize);
        
    for(int i = 0; i < matrixSize; ++i)
    {
        re(i) = generator() * std::sqrt(reEigenvalsRe(i));
        im(i) = generator() * std::sqrt(imEigenvalsRe(i));
    }
    
    Blas_Mat_Vec_Mult(reVecs, re, reRot);
    Blas_Mat_Vec_Mult(imVecs, im, imRot);
    
    alm.Set(lMax, lMax);
    for(int l = 0; l <= lMax; ++l)
    {
        for(int m = 0; m <= l; ++m)
        {
            if(l < lMin)
            {
                alm(l, m) = zero;
                continue;
            }
            
            const int i = index(l, m, lMin);
            check(i < matrixSize, "");
            
            alm(l, m) = xcomplex<double>(reRot(i), imRot(i));
        }
    }
    //output_screen("OK" << std::endl);
    
    if(chi2 != NULL)
    {
        check(dof, "");
        int size1 = index1(lMax, lMax, lMin) + 1;
        LaGenMatDouble wholeMat(size1, size1);
        for(int l1 = lMin; l1 <= lMax; ++l1)
            for(int m1 = -l1; m1 <= l1; ++m1)
                for(int l = lMin; l <= lMax; ++l)
                    for(int m = -l; m <= l; ++m)
                    {
                        int i = index1(l, m, lMin);
                        int j = index1(l1, m1, lMin);
                        check(i < size1, "");
                        check(j < size1, "");
                        wholeMat(i, j) = wholeMatrix.element(l1, m1, l, m);
                    }

        LaVectorLongInt pivot(size1);
        LUFactorizeIP(wholeMat, pivot);
        LaLUInverseIP(wholeMat, pivot);
        
        xcomplex<double> chi2Complex = xcomplex<double>(0, 0);
        for(int l1 = lMin; l1 <= lMax; ++l1)
            for(int m1 = -l1; m1 <= l1; ++m1)
                for(int l = lMin; l <= lMax; ++l)
                    for(int m = -l; m <= l; ++m)
                    {
                        int i = index1(l, m, lMin);
                        int j = index1(l1, m1, lMin);
                        check(i < size1, "");
                        check(j < size1, "");
                        
                        xcomplex<double> almThis, al1m1This;
                        if(m >= 0)
                            almThis = alm(l, m);
                        else
                        {
                            almThis = alm(l, -m).conj();
                            if(m % 2)
                                almThis *= -1.0;
                        }
                        
                        if(m1 >= 0)
                            al1m1This = alm(l1, m1);
                        else
                        {
                            al1m1This = alm(l1, -m1).conj();
                            if(m1 % 2)
                                al1m1This *= -1.0;
                        }
                        chi2Complex += almThis * wholeMat(i, j) * al1m1This.conj();
                    }
        
        check(Math::areEqual(chi2Complex.imag(), 0.0, 1e-5), chi2Complex.imag() << " must be 0");
        (*chi2) = chi2Complex.real();
        (*dof) = size1;
    }
}

void
Simulate::simulateAlm(const std::vector<double>& cl, Alm<xcomplex<double> >& alm, int lMax, time_t seed)
{
    if(seed == 0)
        seed = std::time(0);

    check(!cl.empty(), "");

    if(lMax == 0)
        lMax = cl.size() - 1;

    check(lMax <= cl.size() - 1 || lMax >= 0, "invalid lMax = " << lMax);
    
    const int lMin = 0;
    
    //output_screen("Simulating with seed " << seed << "..." << std::endl)
    boost::variate_generator<boost::mt19937, boost::normal_distribution<> > generator(boost::mt19937(seed), boost::normal_distribution<>(0, 1));
    
    alm.Set(lMax, lMax);
    for(int l = 0; l <= lMax; ++l)
    {
        alm(l, 0) = xcomplex<double>(generator() * std::sqrt(cl[l]), 0.0);
        for(int m = 1; m <= l; ++m)
        {
            const double re = generator() * std::sqrt(cl[l] / 2);
            const double im = generator() * std::sqrt(cl[l] / 2);
            alm(l, m) = xcomplex<double>(re, im);
        }
    }
    //output_screen("OK" << std::endl);
}

void
Simulate::diagonalizeMatrix(const LaGenMatDouble& matrix, LaVectorDouble& eigenvalsRe, LaVectorDouble& eigenvalsIm, LaGenMatDouble& vecs)
{
    const int matrixSize = matrix.rows();
    check(matrix.cols() == matrixSize, "the matrix to be diagonalized is not square");

    eigenvalsRe.resize(matrixSize);
    eigenvalsIm.resize(matrixSize);
    vecs.resize(matrixSize, matrixSize);

    //output_screen("Diagonalizing the covariance matrix..." << std::endl);
    LaEigSolve(matrix, eigenvalsRe, eigenvalsIm, vecs);
    //output_screen("OK" << std::endl);
    
#ifdef CHECKS_ON
    //output_screen("Checking that eigenvalues are real and positive..." << std::endl);
    for(int i = 0; i < matrixSize; ++i)
    {
        check(Math::areEqual(eigenvalsIm(i), 0.0, 1e-20), "The eigenvalue imaginary part " << i << " is supposed to be 0 but it is " << eigenvalsIm(i));
        check(eigenvalsRe(i) >= 0, "The eigenvalue " << i << " is non-positive, it is " << eigenvalsRe(i));
    }
    //output_screen("OK" << std::endl);
    
    //output_screen("Checking that the eigenvectors diagonalize the matrix with correct eigenvalues..." << std::endl);
    LaGenMatDouble mat1(matrixSize, matrixSize), diag(matrixSize, matrixSize);
    Blas_Mat_Mat_Mult(vecs, matrix, mat1, true);
    Blas_Mat_Mat_Mult(mat1, vecs, diag, false);
    
    for(int i = 0; i < matrixSize; ++i)
    {
        for(int j = 0; j < matrixSize; ++j)
        {
            if(i == j)
            {
                check(Math::areEqual(diag(i, i), eigenvalsRe(i), 1e-10), "The diagonal element " << i << " of the matrix " << diag(i, i) << " is not equal to the eigenvalue " << eigenvalsRe(i));
            }
            else
            {
                check(Math::areEqual(diag(i, j), 0.0, 1e-4), "The off diagonal element (" << i << " , " << j << ") of the matrix must be 0 but it is " << diag(i, j));
            }
        }
    }
    //output_screen("OK" << std::endl);
    
    //output_screen("Checking that the eigenvectors are orthogonal..." << std::endl);
    LaGenMatDouble mat2(matrixSize, matrixSize);
    Blas_Mat_Mat_Mult(vecs, vecs, mat2, true);
    
    for(int i = 0; i < matrixSize; ++i)
    {
        for(int j = 0; j < matrixSize; ++j)
        {
            if(i == j)
            {
                check(Math::areEqual(mat2(i, i), 1.0, 1e-10), "The diagonal element " << i << " of the matrix " << mat2(i, i) << " is not equal to 1");
            }
            else
            {
                check(Math::areEqual(mat2(i, j), 0.0, 1e-5), "The off diagonal element (" << i << " , " << j << ") of the matrix must be 0 but it is " << mat2(i, j));
            }
        }
    }
    //output_screen("OK" << std::endl);
#endif
}

void
Simulate::simulateWhiteNoise(Healpix_Map<double>& map, double noiseVal, time_t seed)
{
    if(seed == 0)
        seed = std::time(0);

    boost::variate_generator<boost::mt19937, boost::normal_distribution<> > generator(boost::mt19937(seed), boost::normal_distribution<>(0, noiseVal));
    for(long i = 0; i < map.Npix(); ++i)
        map[i] = generator();
}

