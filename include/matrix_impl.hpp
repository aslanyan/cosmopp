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
Matrix<T>::Matrix(const SymmetricMatrix<DataType>& other)
{
    resize(other.size(), other.size());

#pragma omp parallel for default(shared)
    for(int i = 0; i < rows_; ++i)
        for(int j = 0; j < cols_; ++j)
            (*this)(i, j) = other(i, j);
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
    check(i >= 0 && i < rows(), "invalid index i = " << i << ", should be non-negative and less than " << rows());
    check(j >= 0 && j < cols(), "invalid index j = " << j << ", should be non-negative and less than " << cols());
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
        exceptionStr << "The number of rows and columns must be non-negative! Read " << rows_ << " and " << cols_ << " from " << fileName;
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
Matrix<T>
Matrix<T>::getRow(int i) const
{
    check(i >= 0 && i < rows_, "");
    Matrix<DataType> res(1, cols_);
    for(int j = 0; j < cols_; ++j)
        res(0, j) = (*this)(i, j);
    return res;
}

template<typename T>
Matrix<T>
Matrix<T>::getCol(int i) const
{
    check(i >= 0 && i < cols_, "");
    Matrix<DataType> res(rows_, 1);
    for(int j = 0; j < rows_; ++j)
        res(j, 0) = (*this)(j, i);
    return res;
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
    delete work;

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
int
Matrix<double>::getInverse(Matrix<double>* res) const
{
    res->copy(*this);
    return res->invert();
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

template<typename T>
SymmetricMatrix<T>::SymmetricMatrix(const char* fileName, bool textFile)
{
    if(textFile)
        readFromTextFile(fileName);
    else
        readFromFile(fileName);
}

template<typename T>
const T&
SymmetricMatrix<T>::operator()(int i, int j) const
{
    checkIndices(i, j);
    if(i < j)
        std::swap(i, j);

    return v_[i * (i + 1) / 2 + j];
}

template<typename T>
T&
SymmetricMatrix<T>::operator()(int i, int j)
{
    checkIndices(i, j);
    if(i < j)
        std::swap(i, j);

    return v_[i * (i + 1) / 2 + j];
}

template<typename T>
void
SymmetricMatrix<T>::checkIndices(int i, int j) const
{
    check(i >= 0 && i < size(), "invalid index i = " << i << ", should be non-negative and less than " << size());
    check(j >= 0 && j < size(), "invalid index j = " << j << ", should be non-negative and less than " << size());
}

template<typename T>
void
SymmetricMatrix<T>::resize(int size) 
{
    check(size >= 0, "");

    n_ = size;
    v_.clear();
    v_.resize(n_ * (n_ + 1) / 2);
}

template<typename T>
void
SymmetricMatrix<T>::resize(int size, DataType val) 
{
    check(size >= 0, "");

    n_ = size;
    v_.clear();
    v_.resize(n_ * (n_ + 1) / 2, val);
}

template<typename T>
void
SymmetricMatrix<T>::writeIntoFile(const char* fileName) const
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

    out.write((char*)(&n_), sizeof(n_));
    out.write((char*)(&(v_[0])), v_.size() * sizeof(DataType));
    out.close();
}

template<typename T>
void
SymmetricMatrix<T>::readFromFile(const char* fileName)
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

    in.read((char*)(&n_), sizeof(n_));
    if(n_ < 0)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Invalid matrix size " << n_ << " in the file " << fileName << ". Must be non-negative.";
        exc.set(exceptionStr.str());
        throw exc;
    }

    v_.clear();
    v_.resize(n_ * (n_ + 1) / 2);

    in.read((char*)(&(v_[0])), v_.size() * sizeof(DataType));
    in.close();
}

template<typename T>
void
SymmetricMatrix<T>::writeIntoTextFile(const char* fileName, int precision) const
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

    out << n_ << '\t' << n_ << std::endl;
    if(n_ == 0)
    {
        out.close();
        return;
    }

    check(precision >= 0, "");
    out << std::setprecision(precision);

    for(int i = 0; i < n_; ++i)
    {
        for(int j = 0; j < n_ - 1; ++j)
            out << (*this)(i, j) << '\t';
        out << (*this)(i, n_ - 1) << std::endl;
    }
    out.close();
}

template<typename T>
void
SymmetricMatrix<T>::readFromTextFile(const char* fileName)
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

    int nDummy_;
    in >> n_ >> nDummy_;
    if(n_ < 0)
    {
        std::stringstream exceptionStr;
        exceptionStr << "The matrix size must be non-negative! Read " << n_ << " from " << fileName;
        exc.set(exceptionStr.str());
        throw exc;
    }

    if(nDummy_ != n_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "The matrix in the file " << fileName << " is not symmetric.";
        exc.set(exceptionStr.str());
        throw exc;
    }

    v_.clear();
    v_.resize(n_ * (n_ + 1) / 2);
    for(int i = 0; i < n_; ++i)
    {
        DataType x;
        for(int j = 0; j < i; ++j)
        {
            in >> x;
            if(x != (*this)(j, i))
            {
                std::stringstream exceptionStr;
                exceptionStr << "The matrix in the file is not symmetric! The element (" << i << "," << j << ") is " << x << " while the element (" << j << "," << i << ") is " << (*this)(j, i);
            }
        }
        for(int j = i; j < n_; ++j)
            in >> (*this)(i, j);
    }
    in.close();
}

template<typename T>
void
SymmetricMatrix<T>::copy(const SymmetricMatrix<DataType>& other)
{
    n_ = other.n_;
    v_ = other.v_;
}

template<typename T>
void
SymmetricMatrix<T>::add(const SymmetricMatrix<DataType>& other)
{
    check(n_ == other.n_, "cannot add matrices of different sizes");

#pragma omp parallel for default(shared)
    for(int i = 0; i < v_.size(); ++i)
        v_[i] += other.v_[i];
}

template<typename T>
void
SymmetricMatrix<T>::subtract(const SymmetricMatrix<DataType>& other)
{
    check(n_ == other.n_, "cannot add matrices of different sizes");

#pragma omp parallel for default(shared)
    for(int i = 0; i < v_.size(); ++i)
        v_[i] -= other.v_[i];
}

template<typename T>
void
SymmetricMatrix<T>::multiplyMatrices(const SymmetricMatrix<DataType>& a, const SymmetricMatrix<DataType>& b, SymmetricMatrix<DataType>* res)
{
    check(a.n_ == b.n_, "invalid multiplication, symmetric matrices must have the same size");

    res->resize(a.n_);

#pragma omp parallel for default(shared)
    for(int i = 0; i < a.n_; ++i)
    {
        for(int j = i; j < b.n_; ++j)
        {
            DataType x = 0;
            for(int k = 0; k < a.n_; ++k)
                x += a(i, k) * b(k, j);
            (*res)(i, j) = x;
        }
    }
}

#ifdef COSMO_LAPACK

extern "C"
{
    // Cholesky Factorization
    void dpptrf_(char *uplo, int *n, double *a, int *info);

    // invert from Cholesky factorization
    void dpptri_(char *uplo, int *n, double *a, int *info);

    // tridiagonal reduction
    void dsptrd_(char *uplo, int *n, double *a, double *d, double *e, double *tau, int *info);

    // orthogonal matrix generation after tridiagonal reduction
    void dopgtr_(char *uplo, int *n, double *a, double *tau, double *q, int *ldq, double *work, int *info);

    // eigenvalue and eigenvector calculation
    void dsteqr_(char *compz, int *n, double *d, double *e, double *z, int *ldz, double *work, int *info);

    // eigenvalue and eigenvector calculation positive definite
    void dpteqr_(char *compz, int *n, double *d, double *e, double *z, int *ldz, double *work, int *info);
}

template<>
int
SymmetricMatrix<double>::choleskyFactorize()
{
    check(n_ > 0, "cannot factorize an empty matrix");

    char c = 'U';
    int info;
    int n = n_;

    dpptrf_(&c, &n, &(v_[0]), &info);
    return info;
}

template<>
int
SymmetricMatrix<double>::invertFromCholeskyFactorization()
{
    check(n_ > 0, "matrix is empty");

    char c = 'U';
    int n = n_;
    int info;

    dpptri_(&c, &n, &(v_[0]), &info);
    return info;
}

template<>
int
SymmetricMatrix<double>::invert()
{
    check(n_ > 0, "matrix is empty");
    const int info1 = choleskyFactorize();
    if(info1)
    {
        StandardException exc;
        std::stringstream exceptionStr;
        exceptionStr << "Cholesky factorization failed! info = " << info1 << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    return invertFromCholeskyFactorization();
}

template<>
int
SymmetricMatrix<double>::getInverse(SymmetricMatrix<double>* res) const
{
    res->copy(*this);
    return res->invert();
}

template<>
SymmetricMatrix<double>
SymmetricMatrix<double>::getInverse() const
{
    SymmetricMatrix<double> res;
    getInverse(&res);
    return res;
}

template<>
double
SymmetricMatrix<double>::determinantFromCholeskyFactorization() const
{
    double det = 1;
    for(int i = 0; i < n_; ++i)
        det *= (*this)(i, i);

    det = det * det;

    return det;
}

template<>
double
SymmetricMatrix<double>::determinant() const
{
    check(n_ > 0, "matrix is empty");

    SymmetricMatrix<double> newMat = *this;

    const int info1 = newMat.choleskyFactorize();
    if(info1)
    {
        StandardException exc;
        std::stringstream exceptionStr;
        exceptionStr << "Cholesky factorization failed! info = " << info1 << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    return newMat.determinantFromCholeskyFactorization();
}

template<>
double
SymmetricMatrix<double>::logDetFromCholeskyFactorization(int* sign) const
{
    double logDet = 0;
    *sign = 1;
    for(int i = 0; i < n_; ++i)
    {
        logDet += std::log(std::abs((*this)(i, i)));

        if((*this)(i, i) < 0)
            *sign *= -1;
    }

    return 2 * logDet;
}

template<>
double
SymmetricMatrix<double>::logDet(int* sign) const
{
    check(n_ > 0, "matrix is empty");

    SymmetricMatrix<double> newMat = *this;

    const int info1 = newMat.choleskyFactorize();
    if(info1)
    {
        StandardException exc;
        std::stringstream exceptionStr;
        exceptionStr << "Cholesky factorization failed! info = " << info1 << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    return newMat.logDetFromCholeskyFactorization(sign);
}

template<>
int
SymmetricMatrix<double>::getEigen(std::vector<double>* eigenvals, Matrix<double>* eigenvecs, bool positiveDefinite) const
{
    check(n_ > 0, "matrix is empty");

    std::vector<double> a = v_;
    char uplo = 'U';
    char compz = 'V';
    int n = n_;
    eigenvals->resize(n);
    std::vector<double> e(n - 1);
    std::vector<double> tau(n - 1);
    
    int info;

    dsptrd_(&uplo, &n, &(a[0]), &(eigenvals->at(0)), &(e[0]), &(tau[0]), &info);
    StandardException exc;
    if(info)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Triangular reduction failed! info = " << info << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    eigenvecs->resize(n, n);
    int ldz = n;

    std::vector<double> work(n - 1);

    dopgtr_(&uplo, &n, &(a[0]), &(tau[0]), &(eigenvecs->v_[0]), &ldz, &(work[0]), &info);
    if(info)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Orthogonal matrix generation failed! info = " << info << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    std::vector<double> work1((positiveDefinite ? 4 : 2) * n);
    if(positiveDefinite)
        dpteqr_(&compz, &n, &(eigenvals->at(0)), &(e[0]), &(eigenvecs->v_[0]), &ldz, &(work1[0]), &info);
    else
        dsteqr_(&compz, &n, &(eigenvals->at(0)), &(e[0]), &(eigenvecs->v_[0]), &ldz, &(work1[0]), &info);

    // if positive definite then eigenvals are in descending order, reversing them to still be in ascending order
    if(positiveDefinite)
    {
        for(int i = 0; i < n / 2; ++i)
        {
            std::swap(eigenvals->at(i), eigenvals->at(n - 1 - i));
            for(int j = 0; j < n; ++j)
                std::swap((*eigenvecs)(i, j), (*eigenvecs)(n - 1 - i, j));
        }
    }

    // lapack returns eigenvectors as rows, we want columns
    eigenvecs->transpose();

    return info;
}


#endif

} // namespace Math


#endif

