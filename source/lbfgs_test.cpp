#include <cmath>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <lbfgs.hpp>
#include <cosmo_mpi.hpp>

namespace
{

class LBFGSFunc : public Math::RealFunctionMultiDim
{
public:
    LBFGSFunc(int n) : n_(n)
    {
        check(n_ > 0, "");
    }

    virtual double evaluate(const std::vector<double>& x) const
    {
        check(x.size() == n_, "");

        double res = 0;
        const int processId = CosmoMPI::create().processId();
        for(int i = 0; i < n_; ++i)
            res += std::pow(x[i] - double(i + processId * n_), 2) / (2 * (i + processId * n_ + 1) * (i + processId * n_ + 1));

        double totalRes = 0;
#ifdef COSMO_MPI
        CosmoMPI::create().reduce(&res, &totalRes, 1, CosmoMPI::DOUBLE, CosmoMPI::SUM);
#else
        totalRes = res;
#endif

        return totalRes * totalRes / 2;
    }
private:
    int n_;
};

class LBFGSFuncGrad : public Math::RealFunctionMultiToMulti
{
public:
    LBFGSFuncGrad(int n) : n_(n)
    {
        check(n_ > 0, "");
        sumTag_ = CosmoMPI::create().getCommTag();
    }

    virtual void evaluate(const std::vector<double>& x, std::vector<double>* res) const
    {
        check(x.size() == n_, "");
        res->resize(n_);

        double sum = 0;
        const int processId = CosmoMPI::create().processId();
        for(int i = 0; i < n_; ++i)
            sum += std::pow(x[i] - double(i + processId * n_), 2) / ((i + processId * n_ + 1) * (i + processId * n_ + 1));

        double totalSum = 0;
#ifdef COSMO_MPI
        CosmoMPI::create().reduce(&sum, &totalSum, 1, CosmoMPI::DOUBLE, CosmoMPI::SUM);
#else
        totalSum = sum;
#endif

        if(CosmoMPI::create().isMaster())
        {
            for(int k = 1; k < CosmoMPI::create().numProcesses(); ++k)
                CosmoMPI::create().send(k, &totalSum, 1, CosmoMPI::DOUBLE, sumTag_ + k);
        }
        else
        {
            const int k = CosmoMPI::create().processId();
            check(k != 0, "");
            CosmoMPI::create().recv(0, &totalSum, 1, CosmoMPI::DOUBLE, sumTag_ + k);
        }

        for(int i = 0; i < n_; ++i)
            res->at(i) = (x[i] - double(i + processId * n_)) * totalSum / (2 * (i + processId * n_ + 1) * (i + processId * n_ + 1));
    }

private:
    int n_;
    int sumTag_;
};

void printIter(int iter, double f, double gradNorm, const std::vector<double>& x)
{
    output_screen(iter << '\t' << f << '\t' << gradNorm);
    for(auto it = x.cbegin(); it != x.cend(); ++it)
    {
        output_screen('\t' << *it);
    }
    output_screen(std::endl);
}

} // namespace

int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        using namespace Math;

        output_screen("Input the dimensionality of the problem or it shall be 1 by default." << std::endl);

        int n = 1;
        if(argc > 1)
        {
            std::stringstream str;
            str << argv[1];
            str >> n;

            if(n < 1)
            {
                output_screen("Invalid argument " << argv[1] << " for dimension. Setting it to 1." << std::endl);
                n = 1;
            }
        }

        LBFGSFunc f(n);
        LBFGSFuncGrad g(n);

        LBFGS lbfgs(n, f, g, 10); 

        std::vector<double> x(n, 1000);

        const double epsilon = 1e-3;
        const double gradTol = 1e-3 * n * CosmoMPI::create().numProcesses();

        //lbfgs.minimize(&x, epsilon, 1000000, printIter);
        lbfgs.minimize(&x, epsilon, gradTol, 1000000);

        output_screen("LBFGS is done! Minimum is (found then expected):" << std::endl);
        const int processId = CosmoMPI::create().processId();
        for(int i = 0; i < n; ++i)
        {
            output_screen('\t' << x[i] << '\t' << i + processId * n << std::endl);
        }
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}