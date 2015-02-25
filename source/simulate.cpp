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
#include <random.hpp>
#include <matrix_impl.hpp>
#include <simulate.hpp>

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
    Math::SymmetricMatrix<double> reMatrix(matrixSize), imMatrix(matrixSize);
    
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
    
    std::vector<double> reEigenvalsRe, imEigenvalsRe;
    Math::Matrix<double> reVecs, imVecs;

    reMatrix.getEigen(&reEigenvalsRe, &reVecs, true);
    imMatrix.getEigen(&imEigenvalsRe, &imVecs, true);
    
    //output_screen("Simulating with seed " << seed << "..." << std::endl)
    Math::GaussianGenerator generator(seed, 0, 1);
    
    const xcomplex<double> zero(0, 0);
    Math::Matrix<double> re(matrixSize, 1), im(matrixSize, 1), reRot(matrixSize, 1), imRot(matrixSize, 1);
        
    for(int i = 0; i < matrixSize; ++i)
    {
        re(i, 0) = generator.generate() * std::sqrt(reEigenvalsRe[i]);
        im(i, 0) = generator.generate() * std::sqrt(imEigenvalsRe[i]);
    }
    
    Math::Matrix<double>::multiplyMatrices(reVecs, re, &reRot);
    Math::Matrix<double>::multiplyMatrices(imVecs, im, &imRot);
    
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
            
            alm(l, m) = xcomplex<double>(reRot(i, 0), imRot(i, 0));
        }
    }
    //output_screen("OK" << std::endl);
    
    if(chi2 != NULL)
    {
        check(dof, "");
        int size1 = index1(lMax, lMax, lMin) + 1;
        Math::Matrix<double> wholeMat(size1, size1);
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

        wholeMat.invert();
        
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
    Math::GaussianGenerator generator(seed, 0, 1);
    
    alm.Set(lMax, lMax);
    for(int l = 0; l <= lMax; ++l)
    {
        alm(l, 0) = xcomplex<double>(generator.generate() * std::sqrt(cl[l]), 0.0);
        for(int m = 1; m <= l; ++m)
        {
            const double re = generator.generate() * std::sqrt(cl[l] / 2);
            const double im = generator.generate() * std::sqrt(cl[l] / 2);
            alm(l, m) = xcomplex<double>(re, im);
        }
    }
    //output_screen("OK" << std::endl);
}

void
Simulate::simulateWhiteNoise(Healpix_Map<double>& map, double noiseVal, time_t seed)
{
    if(seed == 0)
        seed = std::time(0);

    Math::GaussianGenerator generator(seed, 0, noiseVal);
    for(long i = 0; i < map.Npix(); ++i)
        map[i] = generator.generate();
}

