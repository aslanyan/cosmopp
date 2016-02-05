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
    int nPar() const;
    void getStarting(std::vector<double> *x) const;
    void getMasses(std::vector<double> *m) const;
    void set(const std::vector<double>& x);
    void get(std::vector<double> *x) const;
    double like() const; // -2ln(like);
    void likeDerivatives(std::vector<double> *d) const; // partial(-2ln(like))/partial(x[i]);
    void output(const std::vector<double>& x, double like);
    bool stop() const;
};
*/

template<typename HMCTraits>
class HMCGeneral
{
public:
    HMCGeneral(HMCTraits *traits, double tauMax, int nMax, int seed = 0);
    ~HMCGeneral();
    void run(int maxIters);

private:
    void generateP();
    double calculatePLike() const;
    double calculateLike() const;

private:
    HMCTraits *traits_;
    int nPar_;
    double tauMax_;
    int nMax_;

    std::unique_ptr<Math::UniformRealGenerator> uniformGen_;
    std::unique_ptr<Math::GaussianGenerator> gaussGen_;

    std::vector<double> mass_;
    std::vector<double> x_, p_, d_;
    std::vector<double> prev_;

    int tauTag_, nTag_, acceptTag_, pLikeTag_, likeTag_;
};

template<typename HMCTraits>
HMCGeneral<HMCTraits>::HMCGeneral(HMCTraits *traits, double tauMax, int nMax, int seed) : traits_(traits), tauMax_(tauMax), nMax_(nMax), nPar_(traits->nPar()), x_(traits->nPar()), p_(traits->nPar()), prev_(traits->nPar()), d_(traits_->nPar())
{
    check(nMax_ >= 1, "");
    check(tauMax_ > 0, "");

    const int gaussSeedTag = CosmoMPI::create().getCommTag();
    tauTag_ = CosmoMPI::create().getCommTag();
    nTag_ = CosmoMPI::create().getCommTag();
    acceptTag_ = CosmoMPI::create().getCommTag();
    pLikeTag_ = CosmoMPI::create().getCommTag();
    likeTag_ = CosmoMPI::create().getCommTag();

    int gaussSeed;

    if(CosmoMPI::create().isMaster())
    {
        int uniformSeed = seed;
        if(uniformSeed == 0)
            uniformSeed = std::time(0);

        gaussSeed = uniformSeed + 1;

        for(int i = 1; i < CosmoMPI::create().numProcesses(); ++i)
        {
            int gaussSeedNew = gaussSeed + i;
            CosmoMPI::create().send(i, &gaussSeedNew, 1, CosmoMPI::INT, gaussSeedTag + i);
        }

        uniformGen_.reset(new Math::UniformRealGenerator(uniformSeed, 0, 1));
        gaussGen_.reset(new Math::GaussianGenerator(gaussSeed, 0, 1));
    }
    else
    {
        const int i = CosmoMPI::create().processId();
        check(i != 0, "");
        CosmoMPI::create().recv(0, &gaussSeed, 1, CosmoMPI::INT, gaussSeedTag + i);
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
void HMCGeneral<HMCTraits>::run(int maxIters)
{
    check(maxIters > 0, "");

    traits_->getStarting(&x_);
    check(x_.size() == nPar_, "");

    traits_->set(x_);
    CosmoMPI::create().barrier();

    double currentLike = calculateLike();

    traits_->output(x_, currentLike);
    
    int total = 0, accepted = 0;
    int gradEvals = 0;

    for(int iter = 0; iter < maxIters; ++iter)
    {
        // x_ should be the current point here, traits_ is set to x_, and currentLike is the likelihood for x_

        double tau;
        int n;
        
        if(CosmoMPI::create().isMaster())
        {
            tau = uniformGen_->generate() * tauMax_;
            n = (int)std::ceil(nMax_ * uniformGen_->generate());
            for(int i = 1; i < CosmoMPI::create().numProcesses(); ++i)
            {
                CosmoMPI::create().send(i, &tau, 1, CosmoMPI::DOUBLE, tauTag_ + i);
                CosmoMPI::create().send(i, &n, 1, CosmoMPI::INT, nTag_ + i);
            }
        }
        else
        {
            const int i = CosmoMPI::create().processId();
            check(i != 0, "");
            CosmoMPI::create().recv(0, &tau, 1, CosmoMPI::DOUBLE, tauTag_ + i);
            CosmoMPI::create().recv(0, &n, 1, CosmoMPI::INT, nTag_ + i);
        }

        check(tau > 0 && tau <= tauMax_, "");
        check(n > 0 && n <= nMax_, "");

        traits_->likeDerivatives(&d_);
        ++gradEvals;
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
            CosmoMPI::create().barrier();
            traits_->likeDerivatives(&d_);
            ++gradEvals;
            for(int k = 0; k < nPar_; ++k)
                p_[k] -= tau / 2 * d_[k] / 2; // divide by 2 because the derivative is of -2ln(like)
        }

        const double newLike = calculateLike();
        const double newPLike = calculatePLike();

        int accept;

        if(CosmoMPI::create().isMaster())
        {
            const double deltaLike = newLike + newPLike - currentLike - oldPLike;
            const double p = std::exp(-deltaLike / 2);
            const double q = uniformGen_->generate();
            accept = (q <= p ? 1 : 0);
            for(int i = 1; i < CosmoMPI::create().numProcesses(); ++i)
            {
                CosmoMPI::create().send(i, &accept, 1, CosmoMPI::INT, acceptTag_ + i);
            }

            ++total;
            if(accept)
                ++accepted;
        }
        else
        {
            const int i = CosmoMPI::create().processId();
            check(i != 0, "");
            CosmoMPI::create().recv(0, &accept, 1, CosmoMPI::INT, acceptTag_ + i);
        }

        CosmoMPI::create().barrier();

        if(accept)
        {
            currentLike = newLike;
        }
        else
        {
            x_.swap(prev_);
            traits_->set(x_);
        }
        CosmoMPI::create().barrier();

        traits_->output(x_, currentLike);

        if(traits_->stop())
            break;
    }

    if(CosmoMPI::create().isMaster())
    {
        output_screen("HMC Sampling finished! Total " << total << " iterations, grad evaluations: " << gradEvals <<", acceptance rate: " << 100 * double(accepted) / total << '%' << std::endl);
    }
}

template<typename HMCTraits>
void HMCGeneral<HMCTraits>::generateP()
{
    check(p_.size() == nPar_, "");
    check(mass_.size() == nPar_, "");

    for(int i = 0; i < nPar_; ++i)
        p_[i] = std::sqrt(mass_[i]) * gaussGen_->generate();
}

template<typename HMCTraits>
double HMCGeneral<HMCTraits>::calculatePLike() const
{
    check(p_.size() == nPar_, "");
    check(mass_.size() == nPar_, "");

    CosmoMPI::create().barrier();

    double res = 0;
    for(int i = 0; i < nPar_; ++i)
        res += p_[i] * p_[i] / mass_[i];

    double totalRes = 0;
#ifdef COSMO_MPI
    CosmoMPI::create().reduce(&res, &totalRes, 1, CosmoMPI::DOUBLE, CosmoMPI::SUM);
#else
    totalRes = res;
#endif

    return totalRes;
}

template<typename HMCTraits>
double HMCGeneral<HMCTraits>::calculateLike() const
{
    double like;
    if(CosmoMPI::create().isMaster())
    {
        like = traits_->like();
        for(int i = 1; i < CosmoMPI::create().numProcesses(); ++i)
        {
            CosmoMPI::create().send(i, &like, 1, CosmoMPI::DOUBLE, likeTag_ + i);
        }
    }
    else
    {
        traits_->like();
        const int i = CosmoMPI::create().processId();
        check(i != 0, "");
        CosmoMPI::create().recv(0, &like, 1, CosmoMPI::DOUBLE, likeTag_ + i);
    }
    return like;
}

} // namespace Math

#endif

