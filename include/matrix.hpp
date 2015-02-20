#ifndef COSMO_PP_MATRIX_HPP
#define COSMO_PP_MATRIX_HPP

#include <vector>

namespace Math
{

template<typename T> class SymmetricMatrix;

template<typename T>
class Matrix
{
    friend class SymmetricMatrix<T>;

public:
    typedef T DataType;

public:
    Matrix() : Matrix(0, 0) {}
    Matrix(int rows, int cols);
    Matrix(int rows, int cols, DataType val);
    Matrix(const char* fileName, bool textFile = false);
    Matrix(const Matrix<DataType>& other) { copy(other); }
    Matrix(const std::vector<DataType>& vec, bool columnVector = true);
    Matrix(const SymmetricMatrix<DataType>& other);

    ~Matrix() {}

    const DataType& operator()(int i, int j) const;
    DataType& operator()(int i, int j);

    int rows() const { return rows_; }
    int cols() const { return cols_; }

    void resize(int rows, int cols);
    void resize(int rows, int cols, DataType val);

    void writeIntoFile(const char* fileName) const;
    void readFromFile(const char* fileName);

    void writeIntoTextFile(const char* fileName, int precision = 3) const;
    void readFromTextFile(const char* fileName);

    Matrix<DataType> getRow(int i) const;
    Matrix<DataType> getCol(int i) const;

    void copy(const Matrix<DataType>& other);

    const Matrix<DataType>& operator=(const Matrix<DataType>& other) { copy(other); return *this; }

    void add(const Matrix<DataType>& other);
    static void addMatrices(const Matrix<DataType>& a, const Matrix<DataType>& b, Matrix<DataType>* res) { res->copy(a); res->add(b); }

    const Matrix<DataType>& operator+=(const Matrix<DataType>& other) { add(other); return *this; }
    Matrix<DataType> operator+(const Matrix<DataType>& other) const { Matrix<DataType> res; addMatrices(*this, other, &res); return res; }

    void subtract(const Matrix<DataType>& other);
    static void subtractMatrices(const Matrix<DataType>& a, const Matrix<DataType>& b, Matrix<DataType>* res) { res->copy(a); res->subtract(b); }

    const Matrix<DataType>& operator-=(const Matrix<DataType>& other) { subtract(other); return *this; }
    Matrix<DataType> operator-(const Matrix<DataType>& other) const { Matrix<DataType> res; subtractMatrices(*this, other, &res); return res; }

    void multiply(const Matrix<DataType>& other) { Matrix<DataType> x; multiplyMatrices(*this, other, &x); copy(x); }
    static void multiplyMatrices(const Matrix<DataType>& a, const Matrix<DataType>& b, Matrix<DataType>* res);

    const Matrix<DataType>& operator*=(const Matrix<DataType>& other) { multiply(other); return *this; }
    Matrix<DataType> operator*(const Matrix<DataType>& other) const { Matrix<DataType> res; multiplyMatrices(*this, other, &res); return res; }

    void getTranspose(Matrix<DataType>* res) const;
    Matrix<DataType> getTranspose() const { Matrix<DataType> res; getTranspose(&res); return res; }
    void transpose() { Matrix<DataType> x; getTranspose(&x); copy(x); }


    int luFactorize(std::vector<int>* pivot);

    int invertFromLUFactorization(std::vector<int> *pivot);
    int invert();

    Matrix<DataType> getInverse() const;
    int getInverse(Matrix<DataType>* res) const;

    double determinantFromLUFactorization(std::vector<int> *pivot) const;
    double determinant() const;

    double logDetFromLUFactorization(std::vector<int> *pivot, int *sign) const;
    double logDet(int *sign) const;

private:
    void checkIndices(int i, int j) const;

private:
    std::vector<DataType> v_;
    int rows_;
    int cols_;
};

template<typename T>
class SymmetricMatrix
{
public:
    typedef T DataType;

public:
    SymmetricMatrix() : SymmetricMatrix(0) {}
    SymmetricMatrix(int size) { resize(size); }
    SymmetricMatrix(int size, DataType val) { resize(size, val); }
    SymmetricMatrix(const char* fileName, bool textFile = false);
    SymmetricMatrix(const SymmetricMatrix<DataType>& other) { copy(other); }

    ~SymmetricMatrix() {}

    const DataType& operator()(int i, int j) const;
    DataType& operator()(int i, int j);

    int size() const { return n_; }

    void resize(int size);
    void resize(int size, DataType val);

    void writeIntoFile(const char* fileName) const;
    void readFromFile(const char* fileName);

    void writeIntoTextFile(const char* fileName, int precision = 3) const;
    void readFromTextFile(const char* fileName);

    void copy(const SymmetricMatrix<DataType>& other);

    const SymmetricMatrix<DataType>& operator=(const SymmetricMatrix<DataType>& other) { copy(other); return *this; }

    void add(const SymmetricMatrix<DataType>& other);
    static void addMatrices(const SymmetricMatrix<DataType>& a, const SymmetricMatrix<DataType>& b, SymmetricMatrix<DataType>* res) { res->copy(a); res->add(b); }

    const SymmetricMatrix<DataType>& operator+=(const SymmetricMatrix<DataType>& other) { add(other); return *this; }
    SymmetricMatrix<DataType> operator+(const SymmetricMatrix<DataType>& other) const { SymmetricMatrix<DataType> res; addMatrices(*this, other, &res); return res; }

    void subtract(const SymmetricMatrix<DataType>& other);
    static void subtractMatrices(const SymmetricMatrix<DataType>& a, const SymmetricMatrix<DataType>& b, SymmetricMatrix<DataType>* res) { res->copy(a); res->subtract(b); }

    const SymmetricMatrix<DataType>& operator-=(const SymmetricMatrix<DataType>& other) { subtract(other); return *this; }
    SymmetricMatrix<DataType> operator-(const SymmetricMatrix<DataType>& other) const { SymmetricMatrix<DataType> res; subtractMatrices(*this, other, &res); return res; }

    void multiply(const SymmetricMatrix<DataType>& other) { SymmetricMatrix<DataType> x; multiplyMatrices(*this, other, &x); copy(x); }
    static void multiplyMatrices(const SymmetricMatrix<DataType>& a, const SymmetricMatrix<DataType>& b, SymmetricMatrix<DataType>* res);

    const SymmetricMatrix<DataType>& operator*=(const SymmetricMatrix<DataType>& other) { multiply(other); return *this; }
    SymmetricMatrix<DataType> operator*(const SymmetricMatrix<DataType>& other) const { SymmetricMatrix<DataType> res; multiplyMatrices(*this, other, &res); return res; }


    int choleskyFactorize();

    int invertFromCholeskyFactorization();
    int invert();

    SymmetricMatrix<DataType> getInverse() const;
    int getInverse(SymmetricMatrix<DataType>* res) const;

    double determinantFromCholeskyFactorization() const;
    double determinant() const;

    double logDetFromCholeskyFactorization(int *sign) const;
    double logDet(int *sign) const;

    int getEigen(std::vector<double>* eigenvals, Matrix<double>* eigenvecs, bool positiveDefinite = false) const;

private:
    void checkIndices(int i, int j) const;

private:
    std::vector<DataType> v_;
    int n_;
};

} // namespace Math

#endif

