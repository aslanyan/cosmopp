#ifndef COSMO_PP_HMC_GENERAL_HPP
#define COSMO_PP_HMC_GENERAL_HPP

#include <ctime>
#include <memory>

#include <macros.hpp>
#include <cosmo_mpi.hpp>
#include <random.hpp>

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
    // multiply with another vector TERM BY TERM
    void multiply(const LargeVector& other);
    // divide by another vector TERM BY TERM
    void divide(const LargeVector& other);
    // take power of elements TERM BY TERM
    void pow(double p);
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

    // generate white noise with given amplitude
    void whitenoise(int seed, LargeVector* x, double amplitude);
};
*/

template <typename LargeVector, typename LargeVectorFactory, typename Function>
class HMCGeneral
{
public:
    HMCGeneral(LargeVectorFactory *factory, Function *f, const LargeVector& starting, const LargeVector& masses, double tauMax, int nMax, int seed = 0);
    ~HMCGeneral();

    // void operator()(const LargeVector& x, double like);
    template<typename CallBack>
    void run(int maxIters, CallBack* callback);

    void stop() { stop_ = true; }

private:
    void generateP();
    double calculatePLike() const;

private:
    LargeVectorFactory *factory_;
    Function *f_;
    double tauMax_;
    int nMax_;
    int gaussSeed_;

    std::unique_ptr<Math::UniformRealGenerator> uniformGen_;

    std::unique_ptr<LargeVector> mass_, massSqrt_;
    std::unique_ptr<LargeVector> x_, p_, d_, prev_, temp_;

    bool stop_;

    CosmoMPI& mpi_;
};

template <typename LargeVector, typename LargeVectorFactory, typename Function>
HMCGeneral<LargeVector, LargeVectorFactory, Function>::HMCGeneral(LargeVectorFactory *factory, Function *f, const LargeVector& starting, const LargeVector& masses, double tauMax, int nMax, int seed) : factory_(factory), f_(f), tauMax_(tauMax), nMax_(nMax), mpi_(CosmoMPI::create()), stop_(false)
{
    check(nMax_ >= 1, "");
    check(tauMax_ > 0, "");

    x_.reset(factory->giveMeOne());
    p_.reset(factory->giveMeOne());
    d_.reset(factory->giveMeOne());
    prev_.reset(factory->giveMeOne());
    mass_.reset(factory->giveMeOne());
    massSqrt_.reset(factory->giveMeOne());
    temp_.reset(factory->giveMeOne());

    if(mpi_.isMaster())
    {
        int uniformSeed = seed;
        if(uniformSeed == 0)
            uniformSeed = std::time(0);

        gaussSeed_ = uniformSeed + 1;

        uniformGen_.reset(new Math::UniformRealGenerator(uniformSeed, 0, 1));
    }
    mpi_.bcast(&gaussSeed_, 1, CosmoMPI::INT);

    mass_->copy(masses);
    x_->copy(starting);
    massSqrt_->copy(*mass_);
    massSqrt_->pow(0.5);
}

template <typename LargeVector, typename LargeVectorFactory, typename Function>
HMCGeneral<LargeVector, LargeVectorFactory, Function>::~HMCGeneral()
{
}

template <typename LargeVector, typename LargeVectorFactory, typename Function>
template<typename CallBack>
void HMCGeneral<LargeVector, LargeVectorFactory, Function>::run(int maxIters, CallBack *callback)
{
    check(maxIters > 0, "");

    f_->set(*x_);
    CosmoMPI::create().barrier();

    double currentLike = f_->value();

    (*callback)(*x_, currentLike);
    
    int total = 0, accepted = 0;
    int gradEvals = 0;

    for(int iter = 0; iter < maxIters; ++iter)
    {
        // x_ should be the current point here, f_ is set to x_, and currentLike is the likelihood for x_

        double tau = -1;
        int n = -1;
        
        if(mpi_.isMaster())
        {
            tau = uniformGen_->generate() * tauMax_;
            n = (int)std::ceil(nMax_ * uniformGen_->generate());
        }
        mpi_.bcast(&tau, 1, CosmoMPI::DOUBLE);
        mpi_.bcast(&n, 1, CosmoMPI::INT);

        check(tau > 0 && tau <= tauMax_, "");
        check(n > 0 && n <= nMax_, "");

        f_->derivative(d_.get());
        ++gradEvals;

        generateP();

        prev_->copy(*x_);
        const double oldPLike = calculatePLike();

        for(int j = 0; j < n; ++j)
        {
            p_->add(*d_, -tau / 2 / 2); // divide by 2 because the derivative is of -2ln(like)
            temp_->copy(*p_);
            temp_->divide(*mass_);
            x_->add(*temp_, tau);

            f_->set(*x_);
            mpi_.barrier();
            f_->derivative(d_.get());
            ++gradEvals;
            p_->add(*d_, -tau / 2 / 2); // divide by 2 because the derivative is of -2ln(like)
        }

        const double newLike = f_->value();
        const double newPLike = calculatePLike();

        int accept;

        if(mpi_.isMaster())
        {
            const double deltaLike = newLike + newPLike - currentLike - oldPLike;
            const double p = std::exp(-deltaLike / 2);
            const double q = uniformGen_->generate();
            accept = (q <= p ? 1 : 0);
            ++total;
            if(accept)
                ++accepted;
        }
        mpi_.bcast(&accept, 1, CosmoMPI::INT);

        CosmoMPI::create().barrier();

        if(accept)
        {
            currentLike = newLike;
        }
        else
        {
            x_->swap(*prev_);
            f_->set(*x_);
        }
        CosmoMPI::create().barrier();

        (*callback)(*x_, currentLike);

        if(stop_)
            break;
    }

    if(mpi_.isMaster())
    {
        output_screen("HMC Sampling finished! Total " << total << " iterations, grad evaluations: " << gradEvals <<", acceptance rate: " << 100 * double(accepted) / total << '%' << std::endl);
    }
}

template <typename LargeVector, typename LargeVectorFactory, typename Function>
void HMCGeneral<LargeVector, LargeVectorFactory, Function>::generateP()
{
    f_->whitenoise(gaussSeed_++, p_.get(), 1.0);
    p_->multiply(*massSqrt_);
}

template <typename LargeVector, typename LargeVectorFactory, typename Function>
double HMCGeneral<LargeVector, LargeVectorFactory, Function>::calculatePLike() const
{
    temp_->copy(*p_);
    temp_->divide(*mass_);
    return temp_->dotProduct(*p_);
}

} // namespace Math

#endif

