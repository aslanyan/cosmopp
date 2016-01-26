#include <cmath>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <lbfgs.hpp>

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
        for(int i = 0; i < n_; ++i)
            res += std::pow(x[i] - double(i), 2) / (2 * (i + 1) * i + 1);

        return res * res / 2;
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
    }

    virtual void evaluate(const std::vector<double>& x, std::vector<double>* res) const
    {
        check(x.size() == n_, "");
        res->resize(n_);

        double sum = 0;
        for(int i = 0; i < n_; ++i)
            sum += std::pow(x[i] - double(i), 2) / ((i + 1) * (i + 1));

        for(int i = 0; i < n_; ++i)
            res->at(i) = (x[i] - double(i)) * sum / (2 * (i + 1) * (i + 1));
    }

private:
    int n_;
};

void printIter(int iter, double gradNorm, const std::vector<double>& x)
{
    output_screen(iter << '\t' << gradNorm);
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

        const double epsilon = 1e-7 / n;

        //lbfgs.minimize(&x, epsilon, 1000000, printIter);
        lbfgs.minimize(&x, epsilon, 1000000);

        /*
        output_screen("LBFGS is done! Minimum is (found then expected):" << std::endl);
        for(int i = 0; i < n; ++i)
        {
            output_screen('\t' << x[i] << '\t' << i << std::endl);
        }
        */
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
