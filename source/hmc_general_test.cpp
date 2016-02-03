#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <memory>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <markov_chain.hpp>
#include <cosmo_mpi.hpp>
#include <hmc_general.hpp>

namespace
{

class TestHMCTraits
{
public:
    TestHMCTraits(int nPar, double mean, double sigma, double mass, double starting, std::string fileRoot) : nPar_(nPar), mean_(mean), sigma_(sigma), mass_(mass), starting_(starting)
    {
        check(nPar_ > 0, "");
        check(mass > 0, "");
        check(sigma > 0, "");

        likeTag_ = CosmoMPI::create().getCommTag();
        outTag_ = CosmoMPI::create().getCommTag();

        std::stringstream fileName;
        fileName << fileRoot;
        if(CosmoMPI::create().numProcesses() > 1)
            fileName << "_" << CosmoMPI::create().processId();

        fileName << ".txt";
        out_.open(fileName.str().c_str());

        StandardException exc;
        if(!out_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Cannot write into file " << fileName.str() << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }
    }

    ~TestHMCTraits()
    {
        out_.close();
    }

    int nPar() const { return nPar_; }
    void getStarting(std::vector<double> *x) const
    {
        x->resize(nPar_);
        for(int i = 0; i < nPar_; ++i)
            x->at(i) = starting_;
    }

    void getMasses(std::vector<double> *m) const
    {
        m->resize(nPar_);
        for(int i = 0; i < nPar_; ++i)
            m->at(i) = mass_;

    }

    void set(const std::vector<double>& x) { x_ = x; }
    void get(std::vector<double> *x) const { *x = x_; }
    double like() const // -2ln(like)
    {
        CosmoMPI::create().barrier();
        double myLike = 0;
        for(int i = 0; i < nPar_; ++i)
        {
            const double delta = x_[i] - mean_;
            myLike += delta * delta / (sigma_ * sigma_);
        }

        double totalLike = 0;
#ifdef COSMO_MPI
        CosmoMPI::create().reduce(&myLike, &totalLike, 1, CosmoMPI::DOUBLE, CosmoMPI::SUM);
#else
        totalLike = myLike;
#endif

        return totalLike;
    }

    void likeDerivatives(std::vector<double> *d) const // partial(-2ln(like))/partial(x[i])
    {
        d->resize(nPar_);
        for(int i = 0; i < nPar_; ++i)
            d->at(i) = 2 * (x_[i] - mean_) / (sigma_ * sigma_);
    }

    void output(const std::vector<double>& x, double like)
    {
        check(x.size() == nPar_, "");
        out_ << "1\t" << like;
        for(int i = 0; i < nPar_; ++i)
            out_ << '\t' << x[i];
        out_ << std::endl;
    }

    bool stop() const
    {
        return false;
    }

private:
    int nPar_;
    double mean_;
    double sigma_;
    double mass_;
    double starting_;
    std::vector<double> x_;

    int likeTag_;
    int outTag_;
    
    std::ofstream out_;
};

}

int main(int argc, char *argv[])
{
    try {
        StandardException exc;

        const int n = 5;
        const int mean = 0;
        const int sigma = 5;
        const double mass = 1;
        const int starting = 1;

        const std::string root = "hmc_general_test";

        TestHMCTraits hmcTraits(n, mean, sigma, mass, starting, root);
        Math::HMCGeneral<TestHMCTraits> hmc(&hmcTraits, 5.0, 10);

        hmc.run(10000);

        const unsigned int thin = 1;
        const unsigned long burnin = 100;

        std::stringstream fileName;
        fileName << root;
        if(CosmoMPI::create().numProcesses() > 1)
            fileName << "_" << CosmoMPI::create().processId();
        fileName << ".txt";

        MarkovChain chain(fileName.str().c_str(), burnin, thin);
        for(int i = 0; i < n; ++i)
        {
            std::unique_ptr<Posterior1D> p(chain.posterior(i, Posterior1D::GAUSSIAN_SMOOTHING));
            const double median = p->median();
            double l, u;
            p->get1SigmaTwoSided(l, u);
            output_screen("Param " << i << ":\t" << median << " + " << u - median << " - " << median - l << std::endl);
        }

    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}

