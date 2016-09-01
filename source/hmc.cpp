#include <sstream>

#include <exception_handler.hpp>
#include <math_constants.hpp>
#include <hmc.hpp>

namespace Math
{

HMC::HMC(int nPar, LikelihoodWithDerivs& like, std::string fileRoot, double tauMax, int nMax, time_t seed) : n_(nPar), like_(like), fileRoot_(fileRoot), tauMax_(tauMax), nMax_(nMax), mass_(nPar, 0), starting_(nPar, 0), derivs_(nPar, 0), currentPoint_(nPar, 0), x_(nPar, 0), p_(nPar, 0), paramNames_(nPar), gen_(seed == 0 ? std::time(0) : seed), uniformDist_(0, 1), gaussDist_(0, 1)
{
    check(n_ > 0, "");
    check(nMax_ >= 1, "");
    check(tauMax_ > 0, "");
}

HMC::~HMC()
{
}

void
HMC::setParam(int i, const std::string& name, double mass, double starting)
{
    check(i >= 0 && i < n_, "invalid index " << i);
    paramNames_[i] = name;
    check(mass > 0, "");
    mass_[i] = mass;
    starting_[i] = starting;
}

void
HMC::run(int iters)
{
    check(iters > 0, "");
    currentPoint_ = starting_;
    currentLike_ = like_.calculate(&(currentPoint_[0]), n_);

    std::stringstream chainFileName;
    chainFileName << fileRoot_ << ".txt";
    outChain_.open(chainFileName.str().c_str());
    StandardException exc;
    if(!outChain_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into the output file " << chainFileName.str() << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    outputCurrent();

    int accepted = 0, total = 0;

    //std::ofstream outXP("xp.txt");

    for(int i = 0; i < iters; ++i)
    {
        const double tau = uniformDist_(gen_) * tauMax_;
        const int n = (int)std::ceil(nMax_ * uniformDist_(gen_));

        output_log1("Iteration " << i << ": tau = " << tau << ", n = " << n << std::endl);

        x_ = currentPoint_;
        calculateDerivs();
        generateP();

        const double oldPLike = calculatePLike();

        /*
        for(int k = 0; k < n_; ++k)
            outXP << x_[k] << '\t';
        for(int k = 0; k < n_ - 1; ++k)
            outXP << p_[k] << '\t';
        outXP << p_[n_ - 1] << std::endl;
        */
        
        check(n > 0, "");
        for(int j = 0; j < n; ++j)
        {
            for(int k = 0; k < n_; ++k)
            {
                p_[k] -= tau / 2 * derivs_[k] / 2; // divide by 2 because the derivative is of -2ln(like)
                x_[k] += tau / mass_[k] * p_[k];
            }
            calculateDerivs();
            for(int k = 0; k < n_; ++k)
                p_[k] -= tau / 2 * derivs_[k] / 2; // divide by 2 because the derivative is of -2ln(like)

            /*
            for(int k = 0; k < n_; ++k)
                outXP << x_[k] << '\t';
            for(int k = 0; k < n_ - 1; ++k)
                outXP << p_[k] << '\t';
            outXP << p_[n_ - 1] << std::endl;
            */
        }

        const double newLike = like_.calculate(&(x_[0]), n_);
        const double newPLike = calculatePLike();

        const double deltaLike = newLike + newPLike - currentLike_ - oldPLike;
        const double p = std::exp(-deltaLike / 2);

        //if(p > 1) p = 1;

        ++total;
        const double q = uniformDist_(gen_);
        if(q <= p)
        {
            output_log1("Iteration " << i << ": New point ACCEPTED" << std::endl);
            ++accepted;
            currentPoint_.swap(x_);
            currentLike_ = newLike;
        }
        else
        {
            output_log1("Iteration " << i << ": New point REJECTED" << std::endl);
        }
        outputCurrent();
    }

    //outXP.close();

    outChain_.close();
    output_screen("Acceptance rate: " << double(accepted) / total * 100 << "%" << std::endl)
}

void
HMC::outputCurrent()
{
    check(outChain_, "");
    outChain_ << 1 << "\t" << currentLike_;
    for(int i = 0; i < n_; ++i)
        outChain_ << "\t" << currentPoint_[i];
    outChain_ << std::endl;
}

void
HMC::generateP()
{
    check(p_.size() == n_, "");
    check(mass_.size() == n_, "");

    for(int i = 0; i < n_; ++i)
        p_[i] = std::sqrt(mass_[i]) * gaussDist_(gen_);
}

double
HMC::calculatePLike()
{
    double res = 0;
    for(int i = 0; i < n_; ++i)
        res += p_[i] * p_[i] / mass_[i];

    return res;
}

void
HMC::calculateDerivs()
{
    for(int i = 0; i < n_; ++i)
        derivs_[i] = like_.calculateDeriv(&(x_[0]), n_, i);
}

} // namespace Math
