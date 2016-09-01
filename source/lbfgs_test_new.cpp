#include <cmath>
#include <fstream>
#include <string>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <lbfgs.hpp>
#include <matrix.hpp>
#include <matrix_impl.hpp>
#include <cosmo_mpi.hpp>

namespace
{

class LBFGSFunc : public Math::RealFunctionMultiDim
{
public:
    LBFGSFunc(double b, const Math::Matrix<double>& a, const Math::Matrix<double>& c) : c_(c), a_(a), b_(b)
    {
    }

    virtual double evaluate(const std::vector<double>& x) const
    {
        Math::Matrix<double> m(x);
        return b_ + (a_ * m + m.getTranspose() * c_ * m)(0, 0);
    }
private:
    const Math::Matrix<double> c_;
    const Math::Matrix<double> a_;
    const double b_;
};

class LBFGSFuncGrad : public Math::RealFunctionMultiToMulti
{
public:
    LBFGSFuncGrad(double b, const Math::Matrix<double>& a, const Math::Matrix<double>& c) : c_(c), a_(a), b_(b)
    {
    }

    virtual void evaluate(const std::vector<double>& x, std::vector<double>* res) const
    {
        Math::Matrix<double> m(x);
        Math::Matrix<double> grad = a_.getTranspose() + c_ * m + c_ * m;
        res->resize(x.size());
        for(int i = 0; i < res->size(); ++i)
            (*res)[i] = grad(i, 0);
    }

private:
    const Math::Matrix<double> c_;
    const Math::Matrix<double> a_;
    const double b_;
};

void callback(int iter, double f, double gradNorm, const std::vector<double>& x, const std::vector<double>& g, const std::vector<double>& z)
{
    check(iter >= 0, "");
    std::stringstream fileNameX, fileNameG, fileNameZ;
    fileNameX << "test_files/lbfgs_test_iter_" << iter << "_x.txt";
    fileNameG << "test_files/lbfgs_test_iter_" << iter << "_g.txt";
    fileNameZ << "test_files/lbfgs_test_iter_" << iter << "_z.txt";
    std::ofstream outX(fileNameX.str().c_str());
    std::ofstream outG(fileNameG.str().c_str());
    std::ofstream outZ(fileNameZ.str().c_str());
    for(int i = 0; i < x.size(); ++i)
    {
        outX << x[i];
        outG << g[i];
        outZ << z[i];
        if(i < x.size() - 1)
        {
            outX << ' ';
            outG << ' ';
            outZ << ' ';
        }
    }
    outX.close();
    outG.close();
    outZ.close();
}

} // namespace

int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        using namespace Math;

        const double b = 5.;
        const int n = 16;
        Math::Matrix<double> a(1, n);
        Math::Matrix<double> c(n, n);
        
        for(int i = 0; i < n; ++i)
        {
            a(0, i) = -double(i);
            for(int j = 0; j < n; ++j)
                c(i, j) = n - std::abs(i - j);
        }

        LBFGSFunc f(b, a, c);
        LBFGSFuncGrad fGrad(b, a, c);

        std::vector<double> starting(n, 2.);

        LBFGS lbfgs(n, f, fGrad, starting, 100);

        std::vector<double> min(n, 0);
        lbfgs.minimize(&min, 1e-20, 1e-100, 1000000, callback);
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
