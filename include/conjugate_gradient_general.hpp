#ifndef COSMO_PP_CONJUGATE_GRADIENT_GENERAL_HPP
#define COSMO_PP_CONJUGATE_GRADIENT_GENERAL_HPP

#include <vector>
#include <memory>

#include <macros.hpp>
#include <line_search.hpp>
#include <cosmo_mpi.hpp>

namespace Math
{

/*
class LargeVector
{
public:
    // copy from other, multiplying with coefficient (for MPI, the correct coefficient should be passed for EVERY process)
    void copy(const LargeVector& other, double c = 1.);
    // set all the elements to 0
    void setToZero();
    // get the norm (for MPI, ALL the processes should get the norm)
    double norm() const;
    // dot product with another vector (for MPI, ALL the processes should get the dot product)
    double dotProduct(const LargeVector& other) const;
    // add another vector with a given coefficient (for MPI, the correct coefficient should be passed for EVERY process)
    void add(const LargeVector& other, double c = 1.);
    // swap
    void swap(LargeVector& other);
};
*/

/*
class LargeVectorFactory
{
public:
    // create a new LargeVector with 0 elements
    // the user is in charge of deleting it
    LargeVector* giveMeOne();
};
*/

/*
class Function
{
public:
    void set(const LargeVector& x);

    // for MPI, ALL the processes should get the function value
    double value();
    void derivative(LargeVector *res);
};
*/

template <typename LargeVector, typename LargeVectorFactory, typename Function>
class CG_General
{
private:
    struct DummyCallBack
    {
        void operator()(int i, double f, double gn, const LargeVector& x, const LargeVector& g, const LargeVector& z)
        {
        }
    };
public:
    enum Method {FLETCHER_REEVES = 0, POLAK_RIBIERE, HESTENES_STIEFEL, DAI_YUAN, METHOD_MAX };
    CG_General(LargeVectorFactory *factory, Function *f, const LargeVector& starting);
    ~CG_General(){}

    void setStarting(const LargeVector& starting);

    // class Callback needs to have
    // void operator()(int iter, double f, double gradNorm, const LargeVector& x, const LargeVector& grad);
    template<typename CallBack>
    double minimize(LargeVector *res, double epsilon = 1e-3, double gNormTol = 1e-5, int maxIter = 1000000, Method m = FLETCHER_REEVES, CallBack* callback = NULL);

    double minimize(LargeVector *res, double epsilon = 1e-3, double gNormTol = 1e-5, int maxIter = 1000000, Method m = FLETCHER_REEVES)
    {
        DummyCallBack* cb = NULL; // this is a hack
        minimize(res, epsilon, gNormTol, maxIter, m, cb);
    }

    void getGradient(LargeVector *g) { g->copy(*g_); }

private:
    Function *f_;
    std::unique_ptr<LargeVector> x_;
    double val_;
    std::unique_ptr<LargeVector> g_, gPrev_;
    std::unique_ptr<LargeVector> z_, zPrev_;
    std::unique_ptr<LargeVector> x0_, g0_;
    double gradNorm_;
    int iter_;

    CosmoMPI& mpi_;
    double rate_;
};

template<typename LargeVector, typename LargeVectorFactory, typename Function>
CG_General<LargeVector, LargeVectorFactory, Function>::CG_General(LargeVectorFactory* factory, Function *f, const LargeVector &starting): f_(f), mpi_(CosmoMPI::create()), rate_(1.0)
{
    // allocate the memory
    x_.reset(factory->giveMeOne());
    g_.reset(factory->giveMeOne());
    gPrev_.reset(factory->giveMeOne());
    z_.reset(factory->giveMeOne());
    zPrev_.reset(factory->giveMeOne());
    x0_.reset(factory->giveMeOne());
    g0_.reset(factory->giveMeOne());

    setStarting(starting);
}

template<typename LargeVector, typename LargeVectorFactory, typename Function>
void
CG_General<LargeVector, LargeVectorFactory, Function>::setStarting(const LargeVector& starting)
{
    x_->copy(starting);

    f_->set(*x_);
    val_ = f_->value();
    f_->derivative(g_.get());
    gradNorm_ = g_->norm();

    gPrev_->copy(*g_);
    zPrev_->setToZero();

    iter_ = 0;
    rate_ = 1.0;
}

template<typename LargeVector, typename LargeVectorFactory, typename Function>
template<typename CallBack>
double
CG_General<LargeVector, LargeVectorFactory, Function>::minimize(LargeVector *res, double epsilon, double gNormTol, int maxIter, Method m, CallBack* callback)
{
    mpi_.barrier();

    check(epsilon > 0, "");
    check(gNormTol >= 0, "");
    check(m >= 0 && m < METHOD_MAX, "");

    int thisIter = 0;
    int functionEval = 0;

    if(callback)
        (*callback)(thisIter, val_, gradNorm_, *x_, *g_, *z_);

    int convergedIters = 0;

    while(true)
    {
        z_->copy(*g_, -1);
        if(iter_ > 0)
        {
            double beta = 0;
            double a, b, c, d;
            switch(m)
            {
            case FLETCHER_REEVES:
                b = gPrev_->dotProduct(*gPrev_);
                check(b > 0, "");
                a = g_->dotProduct(*g_);
                beta = a / b;
                break;
            case POLAK_RIBIERE:
                b = gPrev_->dotProduct(*gPrev_);
                check(b > 0, "");
                a = g_->dotProduct(*g_);
                c = g_->dotProduct(*gPrev_);
                beta = (a - c) / b;
                break;
            case HESTENES_STIEFEL:
                a = g_->dotProduct(*g_);
                b = g_->dotProduct(*gPrev_);
                c = zPrev_->dotProduct(*g_);
                d = zPrev_->dotProduct(*gPrev_);
                check(c - d != 0, "");
                beta = (a - b) / (c - d);
                break;
            case DAI_YUAN:
                a = g_->dotProduct(*g_);
                c = zPrev_->dotProduct(*g_);
                d = zPrev_->dotProduct(*gPrev_);
                check(c - d != 0, "");
                beta = a / (c - d);
                break;
            default:
                check(false, "");
                break;
            }
            if(beta < 0)
            {
                if(mpi_.isMaster())
                {
                    output_screen("CG iteration " << thisIter << ": beta is negative so resetting to steepest descent at this step!" << std::endl);
                }
                beta = 0;
            }

            if(beta > 0)
                z_->add(*zPrev_, beta);
        }

        const double zNorm = z_->norm();
        double zg = -z_->dotProduct(*g_) / zNorm;
        if(zg / gradNorm_ < 0.01)
        {
            if(mpi_.isMaster())
            {
                output_screen("CG iteration " << thisIter << ": Weird stuff! The descent direction does not have a sufficient projection into the gradient. Using steepest descent at this step!" << std::endl);
            }
            z_->copy(*g_, -1);
            zg = 1.0;
        }

        const double oldVal = val_;

        if(iter_ == 0)
            rate_ = 1.0 / gradNorm_;

        // more-thuente line search
        const double ftol = 0.1;
        const double gtol = 0.5;
        const double xtol = 1e-15;
        const double stpmin = 1e-15;
        const double stpmax = 1e15;
        const int maxfev = 100;
        int nfev = 0;
        
        x0_->copy(*x_);

        gPrev_->copy(*g_);
        zPrev_->copy(*z_);


        const int info = moreThuenteSearch(f_, *x0_, val_, *gPrev_, *z_, rate_, ftol, gtol, xtol, stpmin, stpmax, maxfev, x_.get(), g_.get(), nfev);
        check(info != 0, "info needs to be nonzero but it is " << info << ", step = " << rate_ << " iteration: " << iter_);

        functionEval += nfev;
        gradNorm_ = g_->norm();
        
        ++iter_;
        ++thisIter;

        const double deltaVal = std::abs(val_ - oldVal);
        const double valMax = std::max(std::abs(val_), std::abs(oldVal));
        const double ratio = deltaVal / std::max(valMax, 1.0);
        const int minIter = 10;

        if(ratio < epsilon && iter_ >= minIter)
        {
            ++convergedIters;
            if(convergedIters >= 3)
            {
                if(mpi_.isMaster())
                {
                    output_screen("CG has reached the required precision!" << std::endl);
                }
                break;
            }
        }
        else
            convergedIters = 0;

        if(gradNorm_ < gNormTol)
        {
            if(mpi_.isMaster())
            {
                output_screen("CG gradient norm is now below the required tolerance!" << std::endl);
            }
            break;
        }

        if(callback)
            (*callback)(thisIter, val_, gradNorm_, *x_, *g_, *z_);

        if(thisIter > maxIter)
        {
            if(mpi_.isMaster())
            {
                output_screen("CG has reached the maximum number of iterations of " << maxIter << ". Quitting." << std::endl);
            }

            break;
        }
    }

    if(mpi_.isMaster())
    {
        output_screen("CG has converged after " << thisIter << " iterations. Successfully quitting." << std::endl);
        output_screen("Iterations: " << thisIter << ", function evaluations: " << functionEval << ", function value: " << val_ << ", gradient norm: " << gradNorm_ << std::endl);
    }
    res->copy(*x_);
    return val_;
}

} // namespace Math

#endif

