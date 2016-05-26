#ifndef COSMO_PP_LBFGS_HPP
#define COSMO_PP_LBFGS_HPP

#include <memory>

#include <macros.hpp>
#include <function.hpp>
#include <cosmo_mpi.hpp>
#include <lbfgs_general.hpp>

namespace Math
{

/*
class LBFGS
{
public:
    LBFGS(int n, const RealFunctionMultiDim& f, const RealFunctionMultiToMulti& grad, int m = 10);
    ~LBFGS();

    double minimize(std::vector<double> *x, double epsilon = 1e-3, double gNormTol = 1e-5, int maxIter = 1000000, void (*callback)(int, double, double, const std::vector<double>&) = NULL);

private:
    double norm(const std::vector<double>& x) const;

private:
    const RealFunctionMultiDim& f_;
    const RealFunctionMultiToMulti& grad_;
    const int n_;
    const int m_;

    CosmoMPI& mpi_;

    int alphaTag_, betaTag_, z2gTag_, stopTag_, convergedTag_, gradConvergedTag_, ysTag_, yyTag_;
};
*/

class BasicLargeVector
{
public:
    BasicLargeVector(int n) : v_(n, 0) {}

    // copy from other, multiplying with coefficient (for MPI, the correct coefficient should be passed for EVERY process)
    void copy(const BasicLargeVector& other, double c = 1.);
    // set all the elements to 0
    void setToZero();
    // get the norm (for MPI, the master process should get the total norm)
    double norm() const;
    // dot product with another vector (for MPI, the master process should get the total norm)
    double dotProduct(const BasicLargeVector& other) const;
    // add another vector with a given coefficient (for MPI, the correct coefficient should be passed for EVERY process)
    void add(const BasicLargeVector& other, double c = 1.);
    // multiply with another vector TERM BY TERM
    void multiply(const BasicLargeVector& other);
    // divide by another vector TERM BY TERM
    void divide(const BasicLargeVector& other);
    // take power of elements TERM BY TERM
    void pow(double p);
    // swap
    void swap(BasicLargeVector& other);

    std::vector<double>& contents() { return v_;}
    const std::vector<double>& contents() const { return v_;}
    
private:
    std::vector<double> v_;
};

class BasicLargeVectorFactory
{
public:
    BasicLargeVectorFactory(int n) : n_(n)
    {
        check(n_ >= 0, "");
    }

    BasicLargeVector* giveMeOne()
    {
        return new BasicLargeVector(n_);
    }

private:
    int n_;
};

class BasicLBFGSFunc
{
public:
    BasicLBFGSFunc(const RealFunctionMultiDim& f, const RealFunctionMultiToMulti& grad) : f_(f), grad_(grad) {}

    void set(const BasicLargeVector& x) { x_ = x.contents(); }
    double value() { return f_.evaluate(x_); }
    void derivative(BasicLargeVector *res) { grad_.evaluate(x_, &(res->contents())); }

    // generate white noise with given amplitude
    void whitenoise(int seed, BasicLargeVector* x, double amplitude);

private:
    std::vector<double> x_;
    const RealFunctionMultiDim& f_;
    const RealFunctionMultiToMulti& grad_;
};

class BasicLBFGSCallBack
{
public:
    BasicLBFGSCallBack() : cb_(NULL)
    {}

    void set(void (*callback)(int, double, double, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&))
    {
        cb_ = callback;
    }

    void operator()(int iter, double f, double gradNorm, const BasicLargeVector& v, const BasicLargeVector& g, const BasicLargeVector& z)
    {
        if(!cb_)
            return;
        cb_(iter, f, gradNorm, v.contents(), g.contents(), z.contents());
    }
private:
    void (*cb_)(int, double, double, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&);
};

class LBFGS
{
public:
    LBFGS(int n, const RealFunctionMultiDim& f, const RealFunctionMultiToMulti& grad, const std::vector<double>& starting, int m = 10) : factory_(n), f_(f, grad)
    {
        s_ = factory_.giveMeOne();
        s_->contents() = starting;
        lbfgs_.reset(new LBFGS_General<BasicLargeVector, BasicLargeVectorFactory, BasicLBFGSFunc>(&factory_, &f_, *s_, m)); 
    }

    void setStarting(const std::vector<double>& starting)
    {
        s_->contents() = starting;
        lbfgs_->setStarting(*s_);
    }

    double minimize(std::vector<double> *res, double epsilon = 1e-3, double gNormTol = 1e-5, int maxIter = 1000000, void (*callback)(int, double, double, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&) = NULL)
    {
        cb_.set(callback);
        const double val = lbfgs_->minimize(s_, epsilon, gNormTol, maxIter, &cb_);
        (*res) = s_->contents();
        return val;
    }

private:
    BasicLargeVectorFactory factory_;
    BasicLBFGSFunc f_;

    BasicLargeVector *s_;

    std::unique_ptr<LBFGS_General<BasicLargeVector, BasicLargeVectorFactory, BasicLBFGSFunc> > lbfgs_;
    BasicLBFGSCallBack cb_;
};

} // namespace Math

#endif

