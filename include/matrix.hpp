#ifndef COSMO_PP_MATRIX_HPP
#define COSMO_PP_MATRIX_HPP

#include <vector>

#include <macros.hpp>

namespace Math
{

template<typename T>
class Matrix
{
public:
    typedef T DataType;

public:
    Matrix(int rows, int cols);
    Matrix(int rows, int cols, DataType val);

    ~Matrix();

    const DataType& operator()(int i, int j) const;
    DataType& operator()(int i, int j);

    void resize(int rows, int cols);
    void resize(int rows, int cols, DataType val);

private:
    void checkIndices(int i, int j) const;

private:
    std::vector<DataType> v_;
    int rows_;
    int cols_;
};

} // namespace Math

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
Matrix<T>::~Matrix()
{
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

} // namespace Math

#endif

