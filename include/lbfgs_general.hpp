#ifndef COSMO_PP_LBFGS_GENERAL_HPP
#define COSMO_PP_LBFGS_GENERAL_HPP

#include <vector>
#include <memory>

#include <macros.hpp>
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
    // get the norm (for MPI, the master process should get the total norm)
    double norm() const;
    // dot product with another vector (for MPI, the master process should get the total norm)
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
    double value();
    void derivative(LargeVector *res);
}
*/

template <typename LargeVector, typename LargeVectorFactory, typename Function>
class LBFGS_General
{
private:
    struct DummyCallBack
    {
        void operator()(int i, double f, double gn, const LargeVector& x, const LargeVector& g)
        {
        }
    };
public:
    LBFGS_General(LargeVectorFactory *factory, Function *f, const LargeVector& starting, int m = 10);
    ~LBFGS_General(){}

    void setStarting(const LargeVector& starting);

    // class Callback needs to have
    // void operator()(int iter, double f, double gradNorm, const LargeVector& x, const LargeVector& grad);
    template<typename CallBack>
    double minimize(LargeVector *res, double epsilon = 1e-3, double gNormTol = 1e-5, int maxIter = 1000000, CallBack* callback = NULL);

    double minimize(LargeVector *res, double epsilon = 1e-3, double gNormTol = 1e-5, int maxIter = 1000000)
    {
        DummyCallBack* cb = NULL; // this is a hack
        minimize(res, epsilon, gNormTol, maxIter, cb);
    }

    void getGradient(LargeVector *g) { g->copy(*g_); }

private:
    Function *f_;
    std::unique_ptr<LargeVector> x_, xPrev_;
    int m_;
    std::vector<std::unique_ptr<LargeVector> > s_, y_;
    std::vector<double> rho_, alpha_;
    double val_;
    std::unique_ptr<LargeVector> g_, gPrev_;
    std::unique_ptr<LargeVector> q_;
    std::unique_ptr<LargeVector> z_;
    double gradNorm_;
    int iter_;
    std::unique_ptr<LargeVector> searchX_;

    double H0k_;

    CosmoMPI& mpi_;
};

template<typename LargeVector, typename LargeVectorFactory, typename Function>
LBFGS_General<LargeVector, LargeVectorFactory, Function>::LBFGS_General(LargeVectorFactory* factory, Function *f, const LargeVector &starting, int m): f_(f), m_(m), s_(m), y_(m), rho_(m), alpha_(m), mpi_(CosmoMPI::create())
{
    check(m_ > 0, "");
    
    // allocate the memory
    x_.reset(factory->giveMeOne());
    xPrev_.reset(factory->giveMeOne());
    g_.reset(factory->giveMeOne());
    gPrev_.reset(factory->giveMeOne());
    q_.reset(factory->giveMeOne());
    z_.reset(factory->giveMeOne());
    searchX_.reset(factory->giveMeOne());
    for(int i = 0; i < m_; ++i)
    {
        s_[i].reset(factory->giveMeOne());
        y_[i].reset(factory->giveMeOne());
    }

    setStarting(starting);
}

template<typename LargeVector, typename LargeVectorFactory, typename Function>
void
LBFGS_General<LargeVector, LargeVectorFactory, Function>::setStarting(const LargeVector& starting)
{
    x_->copy(starting);

    check(s_.size() == m_, "");
    check(y_.size() == m_, "");
    check(rho_.size() == m_, "");
    check(alpha_.size() == m_, "");
    for(int i = 0; i < m_; ++i)
    {
        s_[i]->setToZero();
        y_[i]->setToZero();
        rho_[i] = 0;
        alpha_[i] = 0;
    }

    f_->set(*x_);
    val_ = f_->value();
    f_->derivative(g_.get());
    gradNorm_ = g_->norm();

    xPrev_->copy(*x_);
    gPrev_->copy(*g_);
    H0k_ = 1;

    iter_ = 0;
}

template<typename LargeVector, typename LargeVectorFactory, typename Function>
template<typename CallBack>
double
LBFGS_General<LargeVector, LargeVectorFactory, Function>::minimize(LargeVector *res, double epsilon, double gNormTol, int maxIter, CallBack* callback)
{
    mpi_.barrier();

    check(epsilon > 0, "");
    check(gNormTol >= 0, "");

    int thisIter = 0;
    int functionEval = 0;

    if(callback)
        (*callback)(thisIter, val_, gradNorm_, *x_, *g_);

    while(true)
    {
        q_->copy(*g_);

        const int m = std::min(m_, iter_); // use this many previous things
        for(int i = 0; i < m; ++i)
        {
            const double dotProduct = s_[i]->dotProduct(*q_);
            alpha_[i] = rho_[i] * dotProduct;
            mpi_.bcast(&(alpha_[i]), 1, CosmoMPI::DOUBLE);
            q_->add(*(y_[i]), -alpha_[i]);
        }
        z_->copy(*q_, H0k_);
        for(int i = m - 1; i >= 0; --i)
        {
            const double dotProduct = y_[i]->dotProduct(*z_);
            double beta = rho_[i] * dotProduct;
            mpi_.bcast(&beta, 1, CosmoMPI::DOUBLE);
            z_->add(*(s_[i]), alpha_[i] - beta);
        }

        double zg = z_->dotProduct(*g_);
        int setZToG = 0;
        if(mpi_.isMaster())
        {
            if(zg <= 0)
            {
                output_screen("LBFGS iteration " << thisIter << ": Weird stuff! The descent direction is not a descent direction. Using conjugate gradient at this step!" << std::endl);
                setZToG = 1;
                zg = gradNorm_ * gradNorm_;
            }
        }
        mpi_.bcast(&setZToG, 1, CosmoMPI::INT);

        if(setZToG)
            z_->copy(*g_);

        const double tau = 0.25, c = 0.01;
        double rate = 1.0;
        searchX_->copy(*x_);
        searchX_->add(*z_, -rate);
        f_->set(*searchX_);
        double newVal = f_->value();
        ++functionEval;
        int searchIter = 0;
        while(true)
        {
            mpi_.barrier();
            int stop = 0;
            if(mpi_.isMaster())
            {
                if(val_ - newVal >= rate * c * zg || searchIter > 100)
                    stop = 1;
            }
            mpi_.bcast(&stop, 1, CosmoMPI::INT);
            if(stop)
                break;

            if(mpi_.isMaster())
            {
                output_screen("Taking a step of size:" << rate * z_->norm() << std::endl);
            }
            rate *= tau;
            searchX_->copy(*x_);
            searchX_->add(*z_, -rate);
            f_->set(*searchX_);
            newVal = f_->value();
            ++functionEval;
            ++searchIter;
        }

        // now move
        x_->copy(*searchX_);
        const double oldVal = val_;
        val_ = newVal;
        f_->derivative(g_.get());
        gradNorm_ = g_->norm();

        const double deltaVal = std::abs(val_ - oldVal);
        const double valMax = std::max(val_, oldVal);
        const double ratio = deltaVal / std::max(valMax, 1.0);
        const int minIter = 10;

        int converged = 0;
        if(mpi_.isMaster())
        {
            if(ratio < epsilon && iter_ >= minIter)
                converged = 1;
        }
        mpi_.bcast(&converged, 1, CosmoMPI::INT);

        if(converged)
        {
            if(mpi_.isMaster())
            {
                output_screen("LBFGS has reached the required precision!" << std::endl);
            }
            break;
        }

        int gradConverged = 0;
        if(mpi_.isMaster())
        {
            if(gradNorm_ < gNormTol)
                gradConverged = 1;
        }
        mpi_.bcast(&gradConverged, 1, CosmoMPI::INT);

        if(gradConverged)
        {
            if(mpi_.isMaster())
            {
                output_screen("LBFGS gradient norm is now below the required tolerance!" << std::endl);
            }
            break;
        }

        // move everything down
        for(int i = m_ - 1; i > 0; --i)
        {
            s_[i]->swap(*(s_[i - 1]));
            y_[i]->swap(*(y_[i - 1]));
            rho_[i] = rho_[i - 1];
        }

        // set the 0 element
        s_[0]->copy(*x_);
        s_[0]->add(*xPrev_, -1);
        y_[0]->copy(*g_);
        y_[0]->add(*gPrev_, -1);
        double ys = s_[0]->dotProduct(*(y_[0]));
        double yy = y_[0]->norm();
        yy = yy * yy;
        mpi_.bcast(&ys, 1, CosmoMPI::DOUBLE);
        mpi_.bcast(&yy, 1, CosmoMPI::DOUBLE);

        if(ys == 0 || yy == 0)
        {
            if(mpi_.isMaster())
            {
                output_screen("LBFGS has not moved for some weird reason! Quitting." << std::endl);
            }

            break;
        }

        rho_[0] = 1 / ys;

        // set H0k
        H0k_ = ys / yy;

        xPrev_->copy(*x_);
        gPrev_->copy(*g_);
        ++iter_;
        ++thisIter;

        if(callback)
            (*callback)(thisIter, val_, gradNorm_, *x_, *g_);

        if(thisIter > maxIter)
        {
            if(mpi_.isMaster())
            {
                output_screen("LBFGS has reached the maximum number of iterations of " << maxIter << ". Quitting." << std::endl);
            }

            break;
        }
    }

    if(mpi_.isMaster())
    {
        output_screen("LBFGS has converged after " << thisIter << " iterations. Successfully quitting." << std::endl);
        output_screen("Iterations: " << thisIter << ", function evaluations: " << functionEval << ", function value: " << val_ << ", gradient norm: " << gradNorm_ << std::endl);
    }
    res->copy(*x_);
    return val_;
}

} // namespace Math

#endif

