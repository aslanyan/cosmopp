#ifndef COSMO_PP_MATRIX_HPP
#define COSMO_PP_MATRIX_HPP

#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>

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

    ~Matrix();

    const DataType& operator()(int i, int j) const;
    DataType& operator()(int i, int j);

    void resize(int rows, int cols);
    void resize(int rows, int cols, DataType val);

    void writeIntoFile(const char* fileName) const;
    void readFromFile(const char* fileName);

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
Matrix<T>::Matrix(const char* fileName, bool textFile)
{
    if(textFile)
    {
    }
    else
        readFromFile(fileName);
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

} // namespace Math

#endif

