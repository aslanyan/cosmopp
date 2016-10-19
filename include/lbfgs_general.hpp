#ifndef COSMO_PP_LBFGS_GENERAL_HPP
#define COSMO_PP_LBFGS_GENERAL_HPP

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

/// A general L-BFGS optimizer.

/// This class can be used to optimize functions (including nonlinear) of possibly very large dimensions using the L-BFGS method.
template <typename LargeVector, typename LargeVectorFactory, typename Function>
class LBFGS_General
{
private:
    struct DummyCallBack
    {
        void operator()(int i, double f, double gn, const LargeVector& x, const LargeVector& g, const LargeVector& z)
        {
        }
    };
public:
    /// Constructor.
    /// \param factory A factory for LargeVector. This needs to be of type that support a function "LargeVector* giveMeOne()" which will create a new LargeVector and return the pointer.
    /// \param f The function to be optimized. This needs to be of type that supports functions "void set(const LargeVector& x)", "double value()", and "void derivative(LargeVector *res)". You first set a point x with set, then you can get the function value and the derivative using value and derivative, respectively.
    /// \param starting The starting point.
    /// \param m The number of previous iterations to store for L-BFGS.
    /// \param moreTheunteLineSearch Determines whether the More-Thuente line search should be used or not. If not, a simple backtracking line search is used.
    LBFGS_General(LargeVectorFactory *factory, Function *f, const LargeVector& starting, int m = 10, bool moreThuenteLineSearch = true);
    
    /// Destructor.
    ~LBFGS_General(){}

    /// Set a new starting point.
    /// This function resets the optimizer completely and sets a new starting point. Can run minimize again after this.
    /// \param starting The new starting point.
    void setStarting(const LargeVector& starting);

    /// Function for minimization (the main function of this class) without callback.
    /// \param res A pointer to LargeVector where the resulting minimum point will be stored.
    /// \param epsilon Threshold for optimization. Will stop if two successive iterations change the function value by less than epsilon.
    /// \param gNormTol Threshold for the norm of the gradient. The minimizer will stop if the gradient norm at the given iteration is less than gNormTol.
    /// \param maxIter Maximum number of iterations. The optimizer will stop if it reaches this maximum number, regardless of convergence.
    /// \return The minimum value found.
    double minimize(LargeVector *res, double epsilon = 1e-3, double gNormTol = 1e-5, int maxIter = 1000000)
    {
        DummyCallBack* cb = NULL; // this is a hack
        return minimize(res, epsilon, gNormTol, maxIter, cb);
    }

    /// Function for minimization (the main function of this class) WITH callback.
    /// \param res A pointer to LargeVector where the resulting minimum point will be stored.
    /// \param epsilon Threshold for optimization. Will stop if two successive iterations change the function value by less than epsilon.
    /// \param gNormTol Threshold for the norm of the gradient. The minimizer will stop if the gradient norm at the given iteration is less than gNormTol.
    /// \param maxIter Maximum number of iterations. The optimizer will stop if it reaches this maximum number, regardless of convergence.
    /// \param callback A pointer to a callback object. This is a template type and needs to support "void operator()(int iter, double f, double gradNorm, const LargeVector& x, const LargeVector& grad)". This operator will be called at each iteration. If callback = NULL then no callback will be used.
    /// \return The minimum value found.
    template<typename CallBack>
    double minimize(LargeVector *res, double epsilon, double gNormTol, int maxIter, CallBack* callback);

    /// Get the gradient at current iteration (can be used after minimization to get the gradient at minimum, for example).
    /// \param g A LargeVector pointer where the result will be returned.
    void getGradient(LargeVector *g) { g->copy(*g_); }

    void replay(LargeVector *x);

    //template<typename LargeVecIterator>
    //void replay(LargeVector *x, LargeVecIterator grad);

    void applyInverseHessian(LargeVector *x);

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
    std::unique_ptr<LargeVector> ss_;
    double gradNorm_;
    int iter_;
    std::unique_ptr<LargeVector> searchX_;

    double H0k_;

    std::vector<double> rates_, H0kSaved_;
    std::vector<bool> usingCGSaved_;

    CosmoMPI& mpi_;
    double rate_;

    const bool moreThuente_;
};

template<typename LargeVector, typename LargeVectorFactory, typename Function>
LBFGS_General<LargeVector, LargeVectorFactory, Function>::LBFGS_General(LargeVectorFactory* factory, Function *f, const LargeVector &starting, int m, bool moreThuenteLineSearch): f_(f), m_(m), s_(m), y_(m), rho_(m), alpha_(m), moreThuente_(moreThuenteLineSearch), mpi_(CosmoMPI::create()), rate_(1.0)
{
    check(m_ > 0, "");
    
    // allocate the memory
    x_.reset(factory->giveMeOne());
    xPrev_.reset(factory->giveMeOne());
    g_.reset(factory->giveMeOne());
    gPrev_.reset(factory->giveMeOne());
    q_.reset(factory->giveMeOne());
    z_.reset(factory->giveMeOne());
    ss_.reset(factory->giveMeOne());
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
    rates_.clear();
    H0kSaved_.clear();
    usingCGSaved_.clear();

    rate_ = 1.0;
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
        (*callback)(thisIter, val_, gradNorm_, *x_, *g_, *z_);

    int convergedIters = 0;

    while(true)
    {
        q_->copy(*g_);

        const int m = std::min(m_, iter_); // use this many previous things
        for(int i = 0; i < m; ++i)
        {
            const double dotProduct = s_[i]->dotProduct(*q_);
            alpha_[i] = rho_[i] * dotProduct;
            q_->add(*(y_[i]), -alpha_[i]);
        }
        H0kSaved_.push_back(H0k_);
        z_->copy(*q_, H0k_);
        for(int i = m - 1; i >= 0; --i)
        {
            const double dotProduct = y_[i]->dotProduct(*z_);
            double beta = rho_[i] * dotProduct;
            z_->add(*(s_[i]), alpha_[i] - beta);
        }

        bool usingCG = false;
        const double zNorm = z_->norm();
        double zg = z_->dotProduct(*g_) / zNorm;
        if(zg / gradNorm_ < 0.01)
        {
            if(mpi_.isMaster())
            {
                output_screen("LBFGS iteration " << thisIter << ": Weird stuff! The descent direction does not have a sufficient projection into the gradient. Using steepest descent at this step!" << std::endl);
            }
            z_->copy(*g_);
            zg = 1.0;
            usingCG = true;
        }

        usingCGSaved_.push_back(usingCG);

        const double oldVal = val_;

        if(iter_ <= 10 || usingCG)
            rate_ = 1.0 / gradNorm_;
        else
            rate_ = 1.0;

        if(moreThuente_)
        {
            const double ftol = 0.1;
            const double gtol = 0.5;
            const double xtol = 1e-15;
            const double stpmin = 1e-15;
            const double stpmax = 1e15;
            const int maxfev = 100;
            int nfev = 0;
            
            ss_->copy(*z_, -1.0);
            searchX_->copy(*x_);
            q_->copy(*g_);


            const int info = moreThuenteSearch(f_, *searchX_, val_, *q_, *ss_, rate_, ftol, gtol, xtol, stpmin, stpmax, maxfev, x_.get(), g_.get(), nfev);
            check(info != 0, "info needs to be nonzero but it is " << info << ", step = " << rate_ << " iteration: " << iter_);

            functionEval += nfev;
            gradNorm_ = g_->norm();
        }
        else
        {
            const double tau = 0.5, c = 1e-5;
            searchX_->copy(*x_);
            searchX_->add(*z_, -rate_);
            f_->set(*searchX_);
            double newVal = f_->value();
            ++functionEval;
            while(true)
            {
                const double valMax = std::max(std::abs(val_), std::abs(newVal));
                if(std::abs(val_ - newVal) / std::max(valMax, 1.0) < epsilon)
                    break;

                if(val_ - newVal >= rate_ * c * zg)
                    break;

                rate_ *= tau;
                searchX_->copy(*x_);
                searchX_->add(*z_, -rate_);
                f_->set(*searchX_);
                newVal = f_->value();
                ++functionEval;
            }

            // now move
            x_->copy(*searchX_);
            val_ = newVal;
            f_->derivative(g_.get());
            gradNorm_ = g_->norm();
        }

        ++iter_;
        ++thisIter;
        rates_.push_back(rate_);

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
                    output_screen("LBFGS has reached the required precision!" << std::endl);
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

        if(callback)
            (*callback)(thisIter, val_, gradNorm_, *x_, *g_, *z_);

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

template<typename LargeVector, typename LargeVectorFactory, typename Function>
void
LBFGS_General<LargeVector, LargeVectorFactory, Function>::replay(LargeVector *x)
{
    check(rates_.size() == iter_, "");
    check(H0kSaved_.size() == iter_, "");
    check(usingCGSaved_.size() == iter_, "");

    check(s_.size() >= iter_ - 1, "don't have enough previous steps saved!");

    for(int it = 0; it < iter_; ++it)
    {
        f_->set(*x);
        const bool usingCG = usingCGSaved_[it];
        if(usingCG)
        {
            f_->derivative(z_.get());
        }
        else
        {
            f_->derivative(q_.get());
            const int m = std::min(m_, it); // use this many previous things
            for(int i = 0; i < m; ++i)
            {
                const double dotProduct = s_[iter_ - 1 - it + i]->dotProduct(*q_);
                alpha_[i] = rho_[iter_ - 1 - it + i] * dotProduct;
                q_->add(*(y_[iter_ - 1 - it + i]), -alpha_[i]);
            }
            H0k_ = H0kSaved_[it];
            z_->copy(*q_, H0k_);
            for(int i = m - 1; i >= 0; --i)
            {
                const double dotProduct = y_[iter_ - 1 - it + i]->dotProduct(*z_);
                double beta = rho_[iter_ - 1 - it + i] * dotProduct;
                z_->add(*(s_[iter_ - 1 - it + i]), alpha_[i] - beta);
            }
        }
        const double rate = rates_[it];
        x->add(*z_, -rate);
    }
}

/*
template<typename LargeVector, typename LargeVectorFactory, typename Function>
template<typename LargeVecIterator>
void
LBFGS_General<LargeVector, LargeVectorFactory, Function>::replay(LargeVector *x, LargeVecIterator grad)
{
    check(rates_.size() == iter_, "");
    check(H0kSaved_.size() == iter_, "");
    check(usingCGSaved_.size() == iter_, "");

    check(s_.size() >= iter_ - 1, "don't have enough previous steps saved!");

    for(int it = 0; it < iter_; ++it)
    {
        const bool usingCG = usingCGSaved_[it];
        if(usingCG)
        {
            z_->copy(*(grad++));
        }
        else
        {
            q_->copy(*(grad++));
            const int m = std::min(m_, it); // use this many previous things
            for(int i = 0; i < m; ++i)
            {
                const double dotProduct = s_[iter_ - 1 - it + i]->dotProduct(*q_);
                alpha_[i] = rho_[iter_ - 1 - it + i] * dotProduct;
                q_->add(*(y_[iter_ - 1 - it + i]), -alpha_[i]);
            }
            H0k_ = H0kSaved_[it];
            z_->copy(*q_, H0k_);
            for(int i = m - 1; i >= 0; --i)
            {
                const double dotProduct = y_[iter_ - 1 - it + i]->dotProduct(*z_);
                double beta = rho_[iter_ - 1 - it + i] * dotProduct;
                z_->add(*(s_[iter_ - 1 - it + i]), alpha_[i] - beta);
            }
        }
        const double rate = rates_[it];
        x->add(*z_, -rate);
    }
}
*/

template<typename LargeVector, typename LargeVectorFactory, typename Function>
void
LBFGS_General<LargeVector, LargeVectorFactory, Function>::applyInverseHessian(LargeVector *x)
{
    q_->copy(*x);
    const int m = std::min(m_, iter_); // use this many previous things
    for(int i = 0; i < m; ++i)
    {
        const double dotProduct = s_[i]->dotProduct(*q_);
        alpha_[i] = rho_[i] * dotProduct;
        q_->add(*(y_[i]), -alpha_[i]);
    }
    x->copy(*q_, H0k_);
    for(int i = m - 1; i >= 0; --i)
    {
        const double dotProduct = y_[i]->dotProduct(*z_);
        double beta = rho_[i] * dotProduct;
        x->add(*(s_[i]), alpha_[i] - beta);
    }
}


} // namespace Math

#endif

