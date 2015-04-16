#ifndef COSMO_PP_MATRIX_HPP
#define COSMO_PP_MATRIX_HPP

#include <macros.hpp>

#include <vector>

namespace Math
{

template<typename T>
class Matrix
{
public:
    typedef T DataType;

public:
    Matrix() : Matrix(0, 0) {}
    Matrix(int rows, int cols);
    Matrix(int rows, int cols, DataType val);
    Matrix(const Matrix<DataType>& other);
    Matrix(const std::vector<DataType>& vec, bool columnVector = true);

    virtual ~Matrix() {}

    virtual const DataType& operator()(int i, int j) const;
    virtual DataType& operator()(int i, int j);

    int rows() const { return rows_; }
    int cols() const { return cols_; }

    virtual void resize(int rows, int cols);
    virtual void resize(int rows, int cols, DataType val);

    virtual void writeIntoFile(const char* fileName) const;
    virtual void readFromFile(const char* fileName);

    virtual void writeIntoTextFile(const char* fileName, int precision = 3) const;
    virtual void readFromTextFile(const char* fileName);

    Matrix<DataType> getRow(int i) const;
    Matrix<DataType> getCol(int i) const;

    virtual void copy(const Matrix<DataType>& other);

    Matrix<DataType>& operator=(const Matrix<DataType>& other) { copy(other); return *this; }

    virtual void add(const Matrix<DataType>& other);
    static void addMatrices(const Matrix<DataType>& a, const Matrix<DataType>& b, Matrix<DataType>* res) { res->copy(a); res->add(b); }

    Matrix<DataType>& operator+=(const Matrix<DataType>& other) { add(other); return *this; }
    Matrix<DataType> operator+(const Matrix<DataType>& other) const { Matrix<DataType> res; addMatrices(*this, other, &res); return res; }

    virtual void subtract(const Matrix<DataType>& other);
    static void subtractMatrices(const Matrix<DataType>& a, const Matrix<DataType>& b, Matrix<DataType>* res) { res->copy(a); res->subtract(b); }

    Matrix<DataType>& operator-=(const Matrix<DataType>& other) { subtract(other); return *this; }
    Matrix<DataType> operator-(const Matrix<DataType>& other) const { Matrix<DataType> res; subtractMatrices(*this, other, &res); return res; }

    void multiply(const Matrix<DataType>& other) { Matrix<DataType> x; multiplyMatrices(*this, other, &x); copy(x); }
    static void multiplyMatrices(const Matrix<DataType>& a, const Matrix<DataType>& b, Matrix<DataType>* res);

    Matrix<DataType>& operator*=(const Matrix<DataType>& other) { multiply(other); return *this; }
    Matrix<DataType> operator*(const Matrix<DataType>& other) const { Matrix<DataType> res; multiplyMatrices(*this, other, &res); return res; }

    virtual void getTranspose(Matrix<DataType>* res) const;
    virtual Matrix<DataType> getTranspose() const { Matrix<DataType> res; getTranspose(&res); return res; }
    virtual void transpose() { Matrix<DataType> x; getTranspose(&x); copy(x); }

    virtual bool isSymmetric() const { return false; }

#ifdef COSMO_LAPACK
    virtual int luFactorize(std::vector<int>* pivot);

    virtual int invertFromLUFactorization(std::vector<int> *pivot);
    virtual int invert();

    int getInverse(Matrix<DataType>* res) const;

    virtual double determinantFromLUFactorization(std::vector<int> *pivot) const;
    virtual double determinant() const;

    virtual double logDetFromLUFactorization(std::vector<int> *pivot, int *sign) const;
    virtual double logDet(int *sign) const;
#endif

protected:
    void checkIndices(int i, int j) const;

protected:
    std::vector<DataType> v_;
    int rows_;
    int cols_;
};

template<typename T>
class SymmetricMatrix : public Matrix<T>
{
public:
    typedef T DataType;
    typedef Matrix<T> BaseType;
    
private:
    using BaseType::rows_;
    using BaseType::cols_;
    using BaseType::v_;
    using BaseType::checkIndices;

public:
    SymmetricMatrix() : BaseType() {}
    SymmetricMatrix(int rows, int cols);
    SymmetricMatrix(int rows, int cols, DataType val);
    SymmetricMatrix(const SymmetricMatrix<DataType>& other);

    virtual ~SymmetricMatrix() {}

    virtual const DataType& operator()(int i, int j) const;
    virtual DataType& operator()(int i, int j);

    virtual void resize(int rows, int cols);
    virtual void resize(int rows, int cols, DataType val);

    virtual void writeIntoFile(const char* fileName) const;
    virtual void readFromFile(const char* fileName);

    virtual void writeIntoTextFile(const char* fileName, int precision = 3) const;
    virtual void readFromTextFile(const char* fileName);

    virtual void copy(const Matrix<DataType>& other);
    virtual void add(const Matrix<DataType>& other);
    virtual void subtract(const Matrix<DataType>& other);

    SymmetricMatrix<DataType> operator+(const SymmetricMatrix<DataType>& other) const { SymmetricMatrix<DataType> res; BaseType::addMatrices(*this, other, &res); return res; }
    SymmetricMatrix<DataType> operator-(const SymmetricMatrix<DataType>& other) const { SymmetricMatrix<DataType> res; subtractMatrices(*this, other, &res); return res; }

    void multiply(const Matrix<DataType>& other) { check(false, "cannot multiply into a symmetric matrix"); }
    Matrix<DataType>& operator*=(const Matrix<DataType>& other) { check(false, "cannot multiply into a symmetric matrix"); return *this; }

    virtual void getTranspose(Matrix<DataType>* res) const { res->copy(*this); }
    virtual Matrix<DataType> getTranspose() const { Matrix<DataType> res; getTranspose(&res); return res; }
    virtual void transpose() { }

    virtual bool isSymmetric() const { return true; }

#ifdef COSMO_LAPACK
    virtual int luFactorize(std::vector<int>* pivot) { check(false, "cannot LU factorize a symmetric matrix"); return -1; }
    virtual int invertFromLUFactorization(std::vector<int> *pivot) { check(false, "cannot LU factorize a symmetric matrix"); return -1; }
    virtual double determinantFromLUFactorization(std::vector<int> *pivot) const { check(false, "cannot LU factorize a symmetric matrix"); return 0; }
    virtual double logDetFromLUFactorization(std::vector<int> *pivot, int *sign) const { check(false, "cannot LU factorize a symmetric matrix"); return 0; }

    int choleskyFactorize();

    int invertFromCholeskyFactorization();
    virtual int invert();

    double determinantFromCholeskyFactorization() const;
    virtual double determinant() const;

    double logDetFromCholeskyFactorization(int *sign) const;
    virtual double logDet(int *sign) const;

    int getEigen(std::vector<double>* eigenvals, Matrix<double>* eigenvecs, bool positiveDefinite = false) const;
#endif
};

} // namespace Math

#endif

