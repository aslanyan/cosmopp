#ifndef COSMO_PP_MATRIX_HPP
#define COSMO_PP_MATRIX_HPP

#include <vector>

namespace Math
{

template<typename T>
class Matrix
{
public:
    typedef T DataType;

public:
    Matrix(int rows = 0, int cols = 0);
    Matrix(int rows, int cols, DataType val);
    Matrix(const char* fileName, bool textFile = false);
    Matrix(const Matrix<DataType>& other) { copy(other); }
    Matrix(const std::vector<DataType>& vec, bool columnVector = true);

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
    void getInverse(Matrix<DataType>* res) const;

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

} // namespace Math

#endif

