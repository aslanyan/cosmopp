#ifndef COSMO_PP_MATRIX_HPP
#define COSMO_PP_MATRIX_HPP

#include <macros.hpp>

#include <vector>

namespace Math
{

/// A general matrix class.
template<typename T>
class Matrix
{
public:
    /// The type of the elements of the matrix.
    typedef T DataType;

public:
    /// Default constructor. Constructs an empty matrix.
    Matrix() : Matrix(0, 0) {}

    /// Constructor. The elements will be initialized to their default values.
    /// \param rows The number of rows.
    /// \param cols The number of columns.
    Matrix(int rows, int cols);

    /// Constructor.
    /// \param rows The number of rows.
    /// \param cols The number of columns.
    /// \param val All of the matrix elements will have this value.
    Matrix(int rows, int cols, DataType val);

    /// Copy constructor.
    /// \param other Another matrix to copy from.
    Matrix(const Matrix<DataType>& other);

    /// Constructor. Creates a single row or column matrix.
    /// \param vec An array of elemnts to initialize the matrix with.
    /// \param columnVector Specifies whether the vector should be a column (true, by default), or a row (false) of the matrix.
    Matrix(const std::vector<DataType>& vec, bool columnVector = true);

    /// Destructor.
    virtual ~Matrix() {}

    /// Element access operator.
    /// \param i The row index.
    /// \param j The column index.
    /// \return Constant reference to the (i, j) element of the matrix.
    virtual const DataType& operator()(int i, int j) const;

    /// Element access operator.
    /// \param i The row index.
    /// \param j The column index.
    /// \return Reference to the (i, j) element of the matrix.
    virtual DataType& operator()(int i, int j);

    /// The number of rows.
    int rows() const { return rows_; }

    /// The number of columns.
    int cols() const { return cols_; }

    /// Resize the matrix. All of the elements will be assigned the default value (any previous values will be erased).
    /// \param rows The new number of rows.
    /// \param cols The new number of columns.
    virtual void resize(int rows, int cols);

    /// Resize the matrix. All of the elements will be assigned the new value val (any previous values will be erased).
    /// \param rows The new number of rows.
    /// \param cols The new number of columns.
    /// \param val All of the elements of the matrix will have this value.
    virtual void resize(int rows, int cols, DataType val);

    /// Write (save) into a binary file.
    /// \param fileName The name of the file.
    virtual void writeIntoFile(const char* fileName) const;

    /// Read (retrieve) from a binary file.
    /// \param fileName The name of the file.
    virtual void readFromFile(const char* fileName);

    /// Write (save) into a text file.
    /// \param fileName The name of the file.
    /// \param precision The precision of the output values.
    virtual void writeIntoTextFile(const char* fileName, int precision = 3) const;

    /// Read (retrieve) from a text file.
    /// \param fileName The name of the file.
    virtual void readFromTextFile(const char* fileName);

    /// Get a given row of the matrix as another matrix.
    /// \param i The index of the row.
    /// \return A new matrix consisting of a single row.
    Matrix<DataType> getRow(int i) const;

    /// Get a given column of the matrix as another matrix.
    /// \param i The index of the column.
    /// \return A new matrix consisting of a single column.
    Matrix<DataType> getCol(int i) const;

    /// Copy from another matrix (any existing values will be erased).
    /// \param other The matrix to copy from.
    virtual void copy(const Matrix<DataType>& other);

    /// Assignment operator.
    /// \param other Another matrix to copy from.
    /// \return A reference to self after the assignment.
    Matrix<DataType>& operator=(const Matrix<DataType>& other) { copy(other); return *this; }

    /// Add another matrix to this matrix (element by element).
    /// \param other The matrix to add.
    virtual void add(const Matrix<DataType>& other);

    /// Add two matrices and write the result into a third one. This is the recommended way to do matrix addition, instead of using the operator +.
    /// \param a The first matrix to add.
    /// \param b The second matrix to add.
    /// \param res A pointer to a matrix where the result will be written.
    static void addMatrices(const Matrix<DataType>& a, const Matrix<DataType>& b, Matrix<DataType>* res) { res->copy(a); res->add(b); }

    /// Add another matrix to the given matrix (element by element).
    /// \param other The other matrix to add.
    /// \return A reference to self after the addition.
    Matrix<DataType>& operator+=(const Matrix<DataType>& other) { add(other); return *this; }

    /// Addition operator. It is recommended to use the addMatrices static function instead since the addition operator returns the result by value which is not efficient.
    /// \param other The right hand side of the operator + (the matrix to add to this).
    /// \return A matrix that's the sum of this and other.
    Matrix<DataType> operator+(const Matrix<DataType>& other) const { Matrix<DataType> res; addMatrices(*this, other, &res); return res; }

    /// Subtract another matrix from the given matrix (element by element).
    /// \param other The other matrix to subtract.
    virtual void subtract(const Matrix<DataType>& other);

    /// Subtract two matrices and write the result into a third one. This is the recommended way to do matrix subtraction, instead of using the operator -.
    /// \param a The matrix to subtract from.
    /// \param b The matrix being subtracted.
    /// \param res A pointer to a matrix where the result will be written.
    static void subtractMatrices(const Matrix<DataType>& a, const Matrix<DataType>& b, Matrix<DataType>* res) { res->copy(a); res->subtract(b); }

    /// Subtract another matrix from the given matrix (element by element).
    /// \param other The other matrix to subtract.
    /// \return A reference to self after the subtraction.
    Matrix<DataType>& operator-=(const Matrix<DataType>& other) { subtract(other); return *this; }

    /// Subtraction operator. It is recommended to use the subtractMatrices static function instead since the subtraction operator returns the result by value which is not efficient.
    /// \param other The other matrix to subtract.
    /// \return A matrix that is the difference of this and other.
    Matrix<DataType> operator-(const Matrix<DataType>& other) const { Matrix<DataType> res; subtractMatrices(*this, other, &res); return res; }

    /// Multiply this matrix with another matrix (this is lhs, other matrix is rhs).
    /// \param other The other matrix.
    void multiply(const Matrix<DataType>& other) { Matrix<DataType> x; multiplyMatrices(*this, other, &x); copy(x); }

    /// Multiply two matrices and write the result into a third one. This is the recommended way to do matrix multiplication, instead of using the operator *.
    /// \param a The left hand side matrix.
    /// \param b The right hand side matrix.
    /// \param res A pointer to a matrix where the result will be written.
    static void multiplyMatrices(const Matrix<DataType>& a, const Matrix<DataType>& b, Matrix<DataType>* res);

    /// Multiply this matrix with another matrix (this is lhs).
    /// \param other The other matrix (rhs).
    /// \return A reference to self after multiplication.
    Matrix<DataType>& operator*=(const Matrix<DataType>& other) { multiply(other); return *this; }

    /// Multiplication operator. It is recommended to use the multiplyMatrices static function instead sincethe multiplication operator returns the result by value which is not efficient.
    /// \param other The other matrix (rhs).
    /// \return A matrix that is the product of this and other.
    Matrix<DataType> operator*(const Matrix<DataType>& other) const { Matrix<DataType> res; multiplyMatrices(*this, other, &res); return res; }

    /// Get the transpose of the matrix.
    /// \param res A pointer to a matrix where the transpose matrix will be returned.
    virtual void getTranspose(Matrix<DataType>* res) const;
    
    /// Get the transpose of the matrix.
    /// \return The transpose matrix.
    virtual Matrix<DataType> getTranspose() const { Matrix<DataType> res; getTranspose(&res); return res; }
    
    /// Transpose the matrix (in place).
    virtual void transpose() { Matrix<DataType> x; getTranspose(&x); copy(x); }

    /// Is this a symmetric matrix. Note that this function does not explicitly check all the elements, it just checks the type (because SymmetricMatrix is a subclass of Matrix). For the Matrix class the result is always false.
    virtual bool isSymmetric() const { return false; }

#ifdef COSMO_LAPACK
    /// LU factorize the matrix (in place).
    /// \param pivot The pivot vector will be returned here (see Lapack documentation).
    /// \return 0 if successful, otherwise an error code (see Lapack documentation).
    virtual int luFactorize(std::vector<int>* pivot);

    /// Invert the matrix. This function should be called after luFactorize.
    /// \param pivot The pivot vector returned by luFactorize.
    /// \return 0 if successful, otherwise an error code (see Lapack documentation).
    virtual int invertFromLUFactorization(std::vector<int> *pivot);

    /// Invert the matrix (in place).
    virtual int invert();

    /// Get the inverse of the matrix.
    /// \param res A pointer to a matrix where the inverse will be returned.
    int getInverse(Matrix<DataType>* res) const;

    /// Calculate the determinant. This function should be called after luFactorize.
    /// \param pivot The pivot vector returned by luFactorize.
    /// \return The determinant of the matrix.
    virtual double determinantFromLUFactorization(std::vector<int> *pivot) const;

    /// Calculate the determinant.
    /// \return The determinant of the matrix.
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

