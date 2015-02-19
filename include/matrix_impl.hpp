#ifndef COSMO_PP_MATRIX_IMPL_HPP
#define COSMO_PP_MATRIX_IMPL_HPP

#ifdef COSMO_OMP
#include <omp.h>
#endif

#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <matrix.hpp>

namespace Math
{

template<typename T>
Matrix<T>::Matrix(int rows, int cols)
{
    resize(rows, cols);
}

template<typename T>
Matrix<T>::Matrix(int rows, int cols, DataType val)
{
    resize(rows, cols, val);
}

template<typename T>
Matrix<T>::Matrix(const char* fileName, bool textFile)
{
    if(textFile)
        readFromTextFile(fileName);
    else
        readFromFile(fileName);
}

template<typename T>
Matrix<T>::Matrix(const std::vector<DataType>& vec, bool columnVector)
{
    v_ = vec;
    if(columnVector)
    {
        rows_ = vec.size();
        cols_ = 1;
    }
    else
    {
        rows_ = 1;
        cols_ = vec.size();
    }
}

template<typename T>
void
Matrix<T>::resize(int rows, int cols)
{
    check(rows >= 0, "");
    check(cols >= 0, "");

    rows_ = rows;
    cols_ = cols;
    v_.clear();
    v_.resize(rows_ * cols_);
}

template<typename T>
void
Matrix<T>::resize(int rows, int cols, DataType val)
{
    check(rows >= 0, "");
    check(cols >= 0, "");

    rows_ = rows;
    cols_ = cols;
    v_.clear();
    v_.resize(rows_ * cols_, val);
}

template<typename T>
const T&
Matrix<T>::operator()(int i, int j) const
{
    checkIndices(i, j);
    return v_[i * cols_ + j];
}

template<typename T>
T&
Matrix<T>::operator()(int i, int j)
{
    checkIndices(i, j);
    return v_[i * cols_ + j];
}

template<typename T>
void
Matrix<T>::checkIndices(int i, int j) const
{
    check(i >= 0 && i < rows_, "invalid index i = " << i << ", should be non-negative and less than " << rows_);
    check(j >= 0 && j < cols_, "invalid index j = " << j << ", should be non-negative and less than " << cols_);
}

template<typename T>
void
Matrix<T>::writeIntoFile(const char* fileName) const
{
    std::ofstream out(fileName, std::ios::binary | std::ios::out);
    StandardException exc;
    if(!out)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into output file " << fileName;
        exc.set(exceptionStr.str());
        throw exc;
    }

    out.write((char*)(&rows_), sizeof(rows_));
    out.write((char*)(&cols_), sizeof(cols_));
    out.write((char*)(&(v_[0])), v_.size() * sizeof(DataType));
    out.close();
}

template<typename T>
void
Matrix<T>::readFromFile(const char* fileName)
{
    std::ifstream in(fileName, std::ios::binary | std::ios::in);
    StandardException exc;
    if(!in)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot read from file " << fileName;
        exc.set(exceptionStr.str());
        throw exc;
    }

    in.read((char*)(&rows_), sizeof(rows_));
    if(rows_ < 0)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Invalid number of rows " << rows_ << " in the file " << fileName << ". Must be non-negative.";
        exc.set(exceptionStr.str());
        throw exc;
    }

    in.read((char*)(&cols_), sizeof(cols_));
    if(cols_ < 0)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Invalid number of columns " << cols_ << " in the file " << fileName << ". Must be non-negative.";
        exc.set(exceptionStr.str());
        throw exc;
    }

    v_.clear();
    v_.resize(rows_ * cols_);

    in.read((char*)(&(v_[0])), v_.size() * sizeof(DataType));
    in.close();
}

template<typename T>
void
Matrix<T>::writeIntoTextFile(const char* fileName, int precision) const
{
    std::ofstream out(fileName);
    StandardException exc;
    if(!out)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into output file " << fileName;
        exc.set(exceptionStr.str());
        throw exc;
    }

    out << rows_ << '\t' << cols_ << std::endl;
    if(rows_ == 0 || cols_ == 0)
    {
        out.close();
        return;
    }

    check(precision >= 0, "");
    out << std::setprecision(precision);

    for(int i = 0; i < rows_; ++i)
    {
        for(int j = 0; j < cols_ - 1; ++j)
            out << (*this)(i, j) << '\t';
        out << (*this)(i, cols_ - 1) << std::endl;
    }
    out.close();
}

template<typename T>
void
Matrix<T>::readFromTextFile(const char* fileName)
{
    std::ifstream in(fileName);
    StandardException exc;
    if(!in)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot read from file " << fileName;
        exc.set(exceptionStr.str());
        throw exc;
    }

    in >> rows_ >> cols_;
    if(rows_ < 0 || cols_ < 0)
    {
        std::stringstream exceptionStr;
        exceptionStr << "The number of rows and columns must be positive! Read " << rows_ << " and " << cols_ << " from " << fileName;
        exc.set(exceptionStr.str());
        throw exc;
    }

    v_.clear();
    v_.resize(rows_ * cols_);
    for(int i = 0; i < rows_; ++i)
    {
        for(int j = 0; j < cols_; ++j)
            in >> (*this)(i, j);
    }
    in.close();
}

template<typename T>
void
Matrix<T>::copy(const Matrix<DataType>& other)
{
    rows_ = other.rows_;
    cols_ = other.cols_;
    v_ = other.v_;
}

template<typename T>
void
Matrix<T>::add(const Matrix<DataType>& other)
{
    check(rows_ == other.rows_, "cannot add matrices of different sizes");
    check(cols_ == other.cols_, "cannot add matrices of different sizes");

#pragma omp parallel for default(shared)
    for(int i = 0; i < v_.size(); ++i)
        v_[i] += other.v_[i];
}

template<typename T>
void
Matrix<T>::subtract(const Matrix<DataType>& other)
{
    check(rows_ == other.rows_, "cannot add matrices of different sizes");
    check(cols_ == other.cols_, "cannot add matrices of different sizes");

#pragma omp parallel for default(shared)
    for(int i = 0; i < v_.size(); ++i)
        v_[i] -= other.v_[i];
}

template<typename T>
void
Matrix<T>::getTranspose(Matrix<DataType>* res) const
{
    res->resize(cols_, rows_);
    for(int i = 0; i < rows_; ++i)
        for(int j = 0; j < cols_; ++j)
            (*res)(j, i) = (*this)(i, j);
}

template<typename T>
void
Matrix<T>::multiplyMatrices(const Matrix<DataType>& a, const Matrix<DataType>& b, Matrix<DataType>* res)
{
    check(a.cols_ == b.rows_, "invalid multiplication, a must have the same number of columns as b rows");

    res->resize(a.rows_, b.cols_);

#pragma omp parallel for default(shared)
    for(int i = 0; i < a.rows_; ++i)
    {
        for(int j = 0; j < b.cols_; ++j)
        {
            DataType x = 0;
            for(int k = 0; k < a.cols_; ++k)
                x += a(i, k) * b(k, j);
            (*res)(i, j) = x;
        }
    }
}

#ifdef COSMO_LAPACK

extern "C"
{
    // LU Factorization
    void dgetrf_(int *m, int *n, double *a, int *lda, int *piv, int *info);

    // invert from LU factorization
    void dgetri_(int *n, double *a, int *lda, int *piv, double *work, int *lwork, int *info);
}

template<>
int
Matrix<double>::luFactorize(std::vector<int>* pivot)
{
    check(rows_ > 0 && cols_ > 0, "cannot factorize an empty matrix");

    pivot->resize(std::min(rows_, cols_));
    int info;
    int m = rows_;
    int n = cols_;
    int lda = rows_;

    dgetrf_(&m, &n, &(v_[0]), &lda, &((*pivot)[0]), &info);
    return info;
}

template<>
int
Matrix<double>::invertFromLUFactorization(std::vector<int>* pivot)
{
    check(rows_ == cols_, "matrix not square");
    check(rows_ > 0, "matrix is empty");

    check(pivot->size() == rows_, "");

    int n = rows_;
    int lda = n;
    int lwork = n * n;
    double *work = new double[lwork];
    int info;

    dgetri_(&n, &(v_[0]), &lda, &((*pivot)[0]), work, &lwork, &info);
    return info;
}

template<>
int
Matrix<double>::invert()
{
    check(rows_ == cols_, "matrix not square");

    std::vector<int> piv;
    const int info1 = luFactorize(&piv);
    if(info1)
    {
        StandardException exc;
        std::stringstream exceptionStr;
        exceptionStr << "LU factorization failed! info = " << info1 << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    return invertFromLUFactorization(&piv);
}

template<>
void
Matrix<double>::getInverse(Matrix<double>* res) const
{
    res->copy(*this);
    const int info = res->invert();
    if(info)
    {
        StandardException exc;
        std::stringstream exceptionStr;
        exceptionStr << "Matrix inversion failed! info = " << info << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
}

template<>
Matrix<double>
Matrix<double>::getInverse() const
{
    Matrix<double> res;
    getInverse(&res);
    return res;
}

template<>
double
Matrix<double>::determinantFromLUFactorization(std::vector<int>* pivot) const
{
    check(rows_ == cols_, "matrix not square");
    check(pivot->size() == rows_, "");

    double det = 1;
    for(int i = 0; i < rows_; ++i)
    {
        det *= (*this)(i, i);
        if((*pivot)[i] != i + 1) // note that pivot elements come from fortran with indices starting from 1
            det *= -1;
    }

    return det;
}

template<>
double
Matrix<double>::determinant() const
{
    check(rows_ == cols_, "matrix not square");

    Matrix<double> newMat = *this;

    std::vector<int> piv;
    const int info1 = newMat.luFactorize(&piv);
    if(info1)
    {
        StandardException exc;
        std::stringstream exceptionStr;
        exceptionStr << "LU factorization failed! info = " << info1 << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    return newMat.determinantFromLUFactorization(&piv);
}

template<>
double
Matrix<double>::logDetFromLUFactorization(std::vector<int>* pivot, int* sign) const
{
    check(rows_ == cols_, "matrix not square");
    check(pivot->size() == rows_, "");

    double logDet = 0;
    *sign = 1;
    for(int i = 0; i < rows_; ++i)
    {
        logDet += std::log(std::abs((*this)(i, i)));

        if((*this)(i, i) < 0)
            *sign *= -1;

        if((*pivot)[i] != i + 1) // note that pivot elements come from fortran with indices starting from 1
            *sign *= -1;
    }

    return logDet;
}

template<>
double
Matrix<double>::logDet(int* sign) const
{
    check(rows_ == cols_, "matrix not square");

    Matrix<double> newMat = *this;

    std::vector<int> piv;
    const int info1 = newMat.luFactorize(&piv);
    if(info1)
    {
        StandardException exc;
        std::stringstream exceptionStr;
        exceptionStr << "LU factorization failed! info = " << info1 << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    return newMat.logDetFromLUFactorization(&piv, sign);
}

#endif

} // namespace Math


#endif

