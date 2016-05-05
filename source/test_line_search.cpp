#include <cstdio>

#include <macros.hpp>
#include <line_search.hpp>
#include <lbfgs.hpp>

using namespace Math;

namespace
{

class MyFunction
{
public:
    MyFunction(const std::vector<double>& x0) : x0_(x0)
    {
        check(!x0_.empty(), "");
    }
    void set(const BasicLargeVector& x)
    {
        x_ = x.contents();
        check(x_.size() == x0_.size(), "");
    }

    // for MPI, ALL the processes should get the function value
    double value()
    {
        const int n = x0_.size();
        check(n > 0, "");
        check(x_.size() == n, "");

        double res = 0;
        for(int i = 0; i < n; ++i)
        {
            const double delta = x_[i] - x0_[i];
            res += delta * delta;
        }
    }

    void derivative(BasicLargeVector *res)
    {
        const int n = x0_.size();
        check(n > 0, "");
        check(x_.size() == n, "");

        std::vector<double> &g = res->contents();
        g.resize(n);
        for(int i = 0; i < n; ++i)
        {
            const double delta = x_[i] - x0_[i];
            g[i] = 2 * delta;
        }
    }

private:
    const std::vector<double> x0_;
    std::vector<double> x_;
};

} // namespace

int main(int argc, char *argv[])
{
    try {
        StandardException exc;

        const int n = 1;
        std::vector<double> x0(n, 0);

        MyFunction func(x0);

        BasicLargeVector p(n);
        for(int i = 0; i < n; ++i)
            p.contents()[i] = 10;

        func.set(p);
        BasicLargeVector g0(n);

        double f = func.value();
        func.derivative(&g0);

        BasicLargeVector s(n);
        s.copy(g0, -2);

        BasicLargeVector x(n), g(n);
        int nFunc = 0;

        double stp = 1;
        const double ftol = 0.5;
        const double gtol = 0.5;
        const double xtol = 0;

        moreThuenteSearch(&func, p, f, g0, s, stp, ftol, gtol, xtol, 0, 10, 100, &x, &g, nFunc);
        std::printf("Step = %.3f\n", stp);
        std::printf("f = %.3f\n", f);
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
