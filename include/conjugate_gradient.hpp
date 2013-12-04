#ifndef COSMO_CPP_CONJUGATE_GRADIENT
#define COSMO_CPP_CONJUGATE_GRADIENT

#include <vector>
#include <cmath>

#include <macros.hpp>

namespace Math
{

/*
struct CGTreats
{
void multiplyByMatrix(const std::vector<double>& original, std::vector<double>& result);
void preconditioner(const std::vector<double>& original, std::vector<double>& result);
};
*/

/// A preconditioned conjugate gradient solver.

/// This class can be used to solve large systems of linear equations using the (preconditioned) conjugate gradient method.
template <class CGTreats>
class ConjugateGradient
{
public:
    /// Constructor.
    /// \param n The number of variables.
    /// \param treats A pointer to an instance of the template treats class. This is used for matrix multiplication and preconditioning. It needs to have a function multiplyByMatrix(const std::vector<double>& original, std::vector<double>& result) for matrix multiplication and preconditioner(const std::vector<double>& original, std::vector<double>& result) for preconditioning. See BasicCGTreats as an example.
    /// \param b The vector b (the system to be solved has the form Ax = b).
    ConjugateGradient(int n, CGTreats* treats, const std::vector<double>& b) : n_(n), treats_(treats), b_(b), r_(n), z_(n), p_(n), x_(n, 0) { check(n_ > 0, ""); check(b_.size() == n_, ""); }
    
    /// The main solver function.
    /// \param precision The precision for solutions.
    /// \param iterations A pointer to an integer where the number of iterations will be written. Ignored if set to NULL (default value).
    /// \return A vector containing the solution.
    const std::vector<double>& solve(double precision = 1e-7, int* iterations = NULL);

private:
    int n_;
    CGTreats* treats_;
    std::vector<double> b_, r_, z_, p_, x_;
};

template <class CGTreats>
const std::vector<double>&
ConjugateGradient<CGTreats>::solve(double precision, int* iterations)
{
    check(precision > 0, "invalid precision = " << precision);

    std::vector<double> ap(n_);

    treats_->multiplyByMatrix(x_, r_);
    for(int i = 0; i < n_; ++i)
        r_[i] = b_[i] - r_[i];

    treats_->preconditioner(r_, z_);
    p_ = z_;

    if(iterations)
        *iterations = 0;
    
    double epsilon;
    do
    {
        double rz = 0;
        for(int i = 0; i < n_; ++i)
            rz += r_[i] * z_[i];

        treats_->multiplyByMatrix(p_, ap);
        double pap = 0;
        for(int i = 0; i < n_; ++i)
            pap += p_[i] * ap[i];

        check(pap != 0, "");

        epsilon = 0;
        const double alpha = rz / pap;
        for(int i = 0; i < n_; ++i)
        {
            x_[i] += alpha * p_[i];
            r_[i] -= alpha * ap[i];

            epsilon += r_[i] * r_[i];
        }

        treats_->preconditioner(r_, z_);
        double rzNew = 0;
        for(int i = 0; i < n_; ++i)
            rzNew += r_[i] * z_[i];

        check(rz != 0, "");
        const double beta = rzNew / rz;
        for(int i = 0; i < n_; ++i)
            p_[i] = z_[i] + beta * p_[i];

        if(iterations)
            ++(*iterations);
    } while(epsilon > precision * precision);

    return x_;
}

/// A simple treats class that can be used as the template parameter of ConjugateGradient.
/// The matrix elements and the preconditioner need to be known and set in advance.
class BasicCGTreats
{
public:
    /// Constructor. Sets the matrix and the preconditioner to a unit matrix.
    /// \param n The number of variables.
    BasicCGTreats(int n) : n_(n), a_(n), p_(n)
    {
        check(n > 0, "");
        for(int i = 0; i < n; ++i)
        {
            a_[i].resize(n, 0);
            p_[i].resize(n, 0);

            a_[i][i] = 1;
            p_[i][i] = 1;
        }
    }

    /// Set the matrix element.
    /// \param i The row (between 0 and n - 1).
    /// \param j The column (between 0 and n - 1).
    /// \param val The value of the element.
    void setMatrix(int i, int j, double val)
    {
        check(i >= 0 && i < n_, "invalid index i = " << i);
        check(j >= 0 && j < n_, "invalid index j = " << j);
        a_[i][j] = val;
    }

    /// Set the matrix element.
    /// \param i The row (between 0 and n - 1).
    /// \param j The column (between 0 and n - 1).
    /// \param val The value of the element.
    void setPreconditioner(int i, int j, double val)
    {
        check(i >= 0 && i < n_, "invalid index i = " << i);
        check(j >= 0 && j < n_, "invalid index j = " << j);
        p_[i][j] = val;
    }

    /// Multiply a given vector by the matrix.
    /// \param x The original vector.
    /// \param y The resulting vector.
    void multiplyByMatrix(const std::vector<double>& x, std::vector<double>& y) const
    {
        check(x.size() == n_, "");
        check(y.size() == n_, "");

        for(int i = 0; i < n_; ++i)
        {
            y[i] = 0;
            for(int j = 0; j < n_; ++j)
                y[i] += a_[i][j] * x[j];
        }
    }

    /// Apply the preconditioner to a given vector.
    /// \param x The original vector.
    /// \param y The resulting vector.
    void preconditioner(const std::vector<double>& x, std::vector<double>& y) const
    {
        check(x.size() == n_, "");
        check(y.size() == n_, "");

        for(int i = 0; i < n_; ++i)
        {
            y[i] = 0;
            for(int j = 0; j < n_; ++j)
                y[i] += p_[i][j] * x[j];
        }
    }

private:
    int n_;
    std::vector<std::vector<double> > a_, p_;
};

} // namespace Math

#endif

