#include <matrix_impl.hpp>

namespace Math
{

#ifdef COSMO_LAPACK

extern "C"
{
    // maltiplication
    int dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *A, int *lda, double *B, int *ldb, double *beta, double *C, int *ldc);

    // LU Factorization
    void dgetrf_(int *m, int *n, double *a, int *lda, int *piv, int *info);

    // invert from LU factorization
    void dgetri_(int *n, double *a, int *lda, int *piv, double *work, int *lwork, int *info);
}

template<>
void
Matrix<double>::multiplyMatrices(const Matrix<double>& a, const Matrix<double>& b, Matrix<double>* res)
{
    check(a.cols_ == b.rows_, "invalid multiplication, a must have the same number of columns as b rows");

    check(!res->isSymmetric(), "the product of two matrices is not necessarily symmetric, even if both are");

    res->resize(a.rows_, b.cols_);

    Matrix<double> &c = *res;

    char transa = 'n';
    char transb = 'n';
    int m = a.rows_;
    int n = b.cols_;
    int k = a.cols_;
    double alpha = 1;
    double beta = 0;
    int lda = a.rows_;
    int ldb = b.rows_;
    int ldc = c.rows_;

    std::vector<double> aVec(a.rows_ * a.cols_), bVec(b.rows_ * b.cols_), cVec(c.rows_ * c.cols_);
    for(int i = 0; i < a.rows_; ++i)
    {
        for(int j = 0; j < a.cols_; ++j)
            aVec[i * a.cols_ + j] = a(i, j);
    }
    for(int i = 0; i < b.rows_; ++i)
    {
        for(int j = 0; j < b.cols_; ++j)
            bVec[i * b.cols_ + j] = b(i, j);
    }

    double *aPt = &(aVec[0]), *bPt = &(bVec[0]);
    const int r = dgemm_(&transa, &transb, &m, &n, &k, &alpha, aPt, &lda, bPt, &ldb, &beta, &(cVec[0]), &ldc);
    check(r == 0, "");

    for(int i = 0; i < c.rows_; ++i)
    {
        for(int j = 0; j < c.cols_; ++j)
            c(i, j) = cVec[i * c.cols_ + j];
    }
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
    check(rows_ == cols_, "");
    check(rows_ > 0, "cannot factorize an empty matrix");

    char c = 'U';
    int info;
    int n = rows_;

    dpptrf_(&c, &n, &(v_[0]), &info);
    return info;
}

template<>
int
SymmetricMatrix<double>::invertFromCholeskyFactorization()
{
    check(rows_ == cols_, "");
    check(rows_ > 0, "matrix is empty");

    char c = 'U';
    int n = rows_;
    int info;

    dpptri_(&c, &n, &(v_[0]), &info);
    return info;
}

template<>
int
SymmetricMatrix<double>::invert()
{
    check(rows_ == cols_, "");
    check(rows_ > 0, "matrix is empty");
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
double
SymmetricMatrix<double>::determinantFromCholeskyFactorization() const
{
    check(rows_ == cols_, "");
    double det = 1;
    for(int i = 0; i < rows_; ++i)
        det *= (*this)(i, i);

    det = det * det;

    return det;
}

template<>
double
SymmetricMatrix<double>::determinant() const
{
    check(rows_ == cols_, "");
    check(rows_ > 0, "matrix is empty");

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
    check(rows_ == cols_, "");
    double logDet = 0;
    *sign = 1;
    for(int i = 0; i < rows_; ++i)
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
    check(rows_ == cols_, "");
    check(rows_ > 0, "matrix is empty");

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
    check(rows_ == cols_, "");
    check(rows_ > 0, "matrix is empty");

    std::vector<double> a = v_;
    char uplo = 'U';
    char compz = 'V';
    int n = rows_;
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

    dopgtr_(&uplo, &n, &(a[0]), &(tau[0]), &((*eigenvecs)(0, 0)), &ldz, &(work[0]), &info);
    if(info)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Orthogonal matrix generation failed! info = " << info << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    std::vector<double> work1((positiveDefinite ? 4 : 2) * n);
    if(positiveDefinite)
        dpteqr_(&compz, &n, &(eigenvals->at(0)), &(e[0]), &((*eigenvecs)(0, 0)), &ldz, &(work1[0]), &info);
    else
        dsteqr_(&compz, &n, &(eigenvals->at(0)), &(e[0]), &((*eigenvecs)(0, 0)), &ldz, &(work1[0]), &info);

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

