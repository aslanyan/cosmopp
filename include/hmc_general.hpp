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
class HMCTraits
{
public:
    unsigned long nPar() const;
    void getStarting(std::vector<double> *x) const;
    void getMasses(std::vector<double> *m) const;
    void set(const std::vector<double>& x);
    void get(std::vector<double> *x) const;
    double like() const; // -2ln(like);
    void likeDerivatives(std::vector<double> *d) const; // partial(-2ln(like))/partial(x[i]);
    void output(const std::vector<double>& x);
};
*/

class HMCGeneralMPIHelper
{
public:
    static void mpiBarrier() const;
    static bool isMaster() const { return processId() == 0; }
    static int processId() const;
    static int numProcesses() const;

    static void sendTo(int i, void *x, int size, int tag) const;
    static void getFrom(int i, void *x, int size, int tag) const;

    static int createNewTag() const;
};

template<typename HMCTraits>
class HMCGeneral
{
public:
    HMCGeneral(HMCTraits *traits, double tauMax, int nMax, time_t seed = 0);
    ~HMCGeneral();
    void run(int iters);

private:
    void mpiBarrier() const;
    bool isMaster() const { return processId() == 0; }
    int processId() const;
    int numProcesses() const;

    void sendTo(int i, void *x, int size, int tag) const;
    void getFrom(int i, void *x, int size, int tag) const;

    int createNewTag() const;

private:
    void generateP();
    double calculatePLike() const;

private:
    HMCTraits *traits_;
    unsigned long nPar_;
    double tauMax_;
    int nMax_;

    std::unique_ptr<Math::UniformRealGenerator> uniformGen_;
    std::unique_ptr<Math::GaussianGenerator> gaussGen_;

    std::vector<double> mass_;
    std::vector<double> x_, p_, d_;
    std::vector<double> prev_;

    int tauTag_, nTag_, acceptTag_, pLikeTag_;
};

template<typename HMCTraits>
HMCGeneral<HMCTraits>::HMCGeneral(HMCTraits *traits, double tauMax, int nMax, time_t seed) : traits_(traits), tauMax_(tauMax), nMax_(nMax), nPar_(traits->nPar()), x_(traits->nPar()), p_(traits->nPar()), prev_(traits->nPar()), d_(traits_->nPar())
{
    check(nMax_ >= 1, "");
    check(tauMax_ > 0, "");

    const int gaussSeedTag = createNewTag();
    tauTag_ = createNewTag();
    nTag_ = createNewTag();
    acceptTag_ = createNewTag();

    time_t gaussSeed;

    if(isMaster())
    {
        time_t uniformSeed = seed;
        if(uniformSeed == 0)
            uniformSeed = std::time(0);

        gaussSeed = uniformSeed + 1;

        for(int i = 1; i < numProcesses(); ++i)
        {
            time_t gaussSeedNew = gaussSeed + i;
            sendTo(i, &gaussSeedNew, sizeof(gaussSeedNew), gaussSeedTag + i);
        }

        uniformGen_.reset(new Math::UniformRealGenerator(uniformSeed, 0, 1));
        gaussGen_.reset(new Math::GaussianGenerator(gaussSeed, 0, 1));
    }
    else
    {
        const int i = processId();
        check(i != 0, "");
        getFrom(0, &gaussSeed, sizeof(gaussSeed), gaussSeedTag + i);
        gaussGen_.reset(new Math::GaussianGenerator(gaussSeed, 0, 1));
    }


    traits_->getMasses(&mass_);
    check(mass_.size() == nPar_, "");
}

template<typename HMCTraits>
HMCGeneral<HMCTraits>::~HMCGeneral()
{
}

template<typename HMCTraits>
void HMCGeneral<HMCTraits>::run(int iters)
{
    check(iters > 0, "");

    traits_->getStarting(&x_);
    check(x_.size() == nPar_, "");

    traits_->set(x_);
    mpiBarrier();
    double currentLike = traits_->like();

    traits_->output(x_);

    for(int iter = 0; iter < iters; ++iter)
    {
        // x_ should be the current point here, traits_ is set to x_, and currentLike is the likelihood for x_

        double tau;
        int n;
        
        if(isMaster())
        {
            tau = uniformGen_->generate(); * tauMax_;
            n = (int)std::ceil(nMax_ * uniformGen_->generate());
            for(int i = 1; i < numProcesses(); ++i)
            {
                sendTo(i, &tau, sizeof(tau), tauTag_ + i);
                sendTo(i, &n, sizeof(n), nTag_ + i);
            }
        }
        else
        {
            const int i = processId();
            check(i != 0, "");
            getFrom(0, &tau, sizeof(tau), tauTag_ + i);
            getFrom(0, &n, sizeof(n), nTag_ + i);
        }

        check(tau > 0 && tau <= tauMax_, "");
        check(n > 0 && n <= nMax_, "");

        traits_->likeDerivatives(&d_);
        check(d_.size() == nPar_, "");

        generateP();

        prev_ = x_;
        const double oldPLike = calculatePLike();

        for(int j = 0; j < n; ++j)
        {
            for(int k = 0; k < nPar_; ++k)
            {
                p_[k] -= tau / 2 * d_[k] / 2; // divide by 2 because the derivative is of -2ln(like)
                x_[k] += tau / mass_[k] * p_[k];
            }
            traits_->set(x_);
            mpiBarrier();
            traits_->likeDerivatives(&d_);
            for(int k = 0; k < n_; ++k)
                p_[k] -= tau / 2 * d_[k] / 2; // divide by 2 because the derivative is of -2ln(like)
        }

        const double newLike = traits_->like();
        const double newPLike = calculatePLike();

        bool accept;

        if(isMaster())
        {
            const double deltaLike = newLike + newPLike - currentLike_ - oldPLike;
            const double p = std::exp(-deltaLike / 2);
            const double q = uniformGen_->generate();
            accept = (q <= p);
            for(int i = 1; i < numProcesses(); ++i)
            {
                sendTo(i, &accept, sizeof(accept), acceptTag_ + i);
            }
        }
        else
        {
            const int i = processId();
            check(i != 0, "");
            getFrom(0, &accept, sizeof(accept), acceptTag_ + i);
        }

        mpiBarrier();

        if(accept)
        {
            currentLike_ = newLike;
        }
        else
        {
            x_.swap(prev_);
            traits_->set(x_);
        }
        mpiBarrier();

        traits_->output(x_);
    }
}

template<typename HMCTraits>
void HMCGeneral<HMCTraits>::generateP()
{
    check(p_.size() == nPar_, "");
    check(mass_.size() == nPar_, "");

    for(int i = 0; i < n_; ++i)
        p_[i] = std::sqrt(mass_[i]) * gaussGen_->generate();
}

template<typename HMCTraits>
double HMCGeneral<HMCTraits>::calculatePLike()
{
    check(p_.size() == nPar_, "");
    check(mass_.size() == nPar_, "");

    mpiBarrier();

    double res = 0;
    for(int i = 0; i < nPar_; ++i)
        res += p_[i] * p_[i] / mass_[i];

    if(isMaster())
    {
        for(int i = 1; i < numProcesses(); ++i)
        {
            double otherRes;
            getFrom(i, &otherRes, sizeof(otherRes), pLikeTag_ + i);
            res += otherRes;
        }
    }
    else
    {
        const int i = processId();
        check(i != 0, "");
        sendTo(0, &res, sizeof(res), pLikeTag_ + i);
        return 0;
    }
}

} // namespace Math

#endif

