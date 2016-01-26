#include <cmath>
#include <algorithm>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <lbfgs.hpp>

#ifdef COSMO_OMP
#include <omp.h>
#endif

namespace Math
{

LBFGS::LBFGS(int n, const RealFunctionMultiDim& f, const RealFunctionMultiToMulti& grad, int m) : n_(n), f_(f), grad_(grad), m_(m)
{
    check(n_ > 0, "the function should have a positive number of parameters");
    check(m_ > 0, "");
}

LBFGS::~LBFGS()
{
}


double
LBFGS::minimize(std::vector<double> *x, double epsilon, int maxIter, void (*callback)(int, double, const std::vector<double>&))
{
    check(x->size() == n_, "");
    check(epsilon > 0, "");

    std::vector<std::vector<double> > s(m_), y(m_);
    for(int i = 0; i < m_; ++i)
    {
        s[i].resize(n_, 0);
        y[i].resize(n_, 0);
    }

    std::vector<double> rho(m_, 0), alpha(m_, 0);

    int iter = 0;

    double val = f_.evaluate(*x);
    std::vector<double> g;
    grad_.evaluate(*x, &g);
    check(g.size() == n_, "");
    double gradNorm = norm(g);

    std::vector<double> xPrev(*x), gPrev(g);
    std::vector<double> q(n_);
    std::vector<double> z(n_);

    std::vector<double> searchX(n_);

    double H0k = 1;
    if(callback)
        callback(iter, gradNorm, *x);

    while(gradNorm > epsilon)
    {
        if(iter > maxIter)
        {
            output_screen("LBFGS has reached the maximum number of iterations of " << maxIter << ". Quitting." << std::endl);
            return gradNorm;
        }

#pragma omp parallel for default(shared)
        for(int i = 0; i < n_; ++i)
            q[i] = g[i];

        const int m = std::min(m_, iter); // use this many previous things
        for(int i = 0; i < m; ++i)
        {
            int j;
            double dotProduct = 0;
#pragma omp parallel for default(shared) private(j) reduction(+:dotProduct)
            for(j = 0; j < n_; ++j)
                dotProduct += s[i][j] * q[j];
            alpha[i] = rho[i] * dotProduct;

#pragma omp parallel for default(shared) private(j)
            for(j = 0; j < n_; ++j)
                q[j] -= alpha[i] * y[i][j];
        }

#pragma omp parallel for default(shared)
        for(int i = 0; i < n_; ++i)
            z[i] = H0k * q[i];

        for(int i = m - 1; i >= 0; --i)
        {
            int j;
            double dotProduct = 0;
#pragma omp parallel for default(shared) private(j) reduction(+:dotProduct)
            for(j = 0; j < n_; ++j)
                dotProduct += y[i][j] * z[j];
            const double beta = rho[i] * dotProduct;

#pragma omp parallel for default(shared) private(j)
            for(j = 0; j < n_; ++j)
                z[j] += (alpha[i] - beta) * s[i][j];
        }

        double rate = 1.0;
        double zg = 0;
        int i;
#pragma omp parallel for default(shared) private(i) reduction(+:zg)
        for(i = 0; i < n_; ++i)
            zg += z[i] * g[i];

        if(zg <= 0)
        {
            output_screen("LBFGS iteration " << iter << ": Weird stuff! The descent direction is not a descent direction. Using conjugate gradient at this step!" << std::endl);
            z = g;
            zg = gradNorm * gradNorm;
        }

        const double tau = 0.5, c = 0.1;

#pragma omp parallel for default(shared) private(i)
        for(i = 0; i < n_; ++i)
            searchX[i] = x->at(i) - rate * z[i];
        double newVal = f_.evaluate(searchX);
        int searchIter = 0;
        while(val - newVal < rate * c * zg)
        {
            if(searchIter > 1000)
                break;
            rate *= tau;

#pragma omp parallel for default(shared)
            for(int i = 0; i < n_; ++i)
                searchX[i] = x->at(i) - rate * z[i];
            newVal = f_.evaluate(searchX);
            ++searchIter;
        }

        // now move
#pragma omp parallel for default(shared) private(i)
        for(i = 0; i < n_; ++i)
            (*x)[i] -= rate * z[i];
        grad_.evaluate(*x, &g);
        check(g.size() == n_, "");
        gradNorm = norm(g);

        // move everything down
        for(i = m_ - 1; i > 0; --i)
        {
            std::swap(s[i], s[i - 1]);
            std::swap(y[i], y[i - 1]);
            rho[i] = rho[i - 1];
        }

        // set the 0 element
        double ys = 0, yy = 0;
#pragma omp parallel for default(shared) private(i)
        for(i = 0; i < n_; ++i)
        {
            s[0][i] = x->at(i) - xPrev[i];
            y[0][i] = g[i] - gPrev[i];
        }
#pragma omp parallel for default(shared) private(i) reduction(+:ys)
        for(i = 0; i < n_; ++i)
        {
            ys += s[0][i] * y[0][i];
        }
#pragma omp parallel for default(shared) private(i) reduction(+:yy)
        for(i = 0; i < n_; ++i)
        {
            yy += y[0][i] * y[0][i];
        }
        if(ys == 0 || yy == 0)
        {
            output_screen("LBFGS has not moved for some weird reason! Quitting." << std::endl);
            return gradNorm;
        }
        rho[0] = 1 / ys;

        // set H0k
        H0k = ys / yy;

#pragma omp parallel for default(shared) private(i)
        for(i = 0; i < n_; ++i)
        {
            xPrev[i] = x->at(i);
            gPrev[i] = g[i];
        }
        ++iter;

        if(callback)
            callback(iter, gradNorm, *x);
    }

    output_screen("LBFGS has converged after " << iter << " iterations. The norm of the gradient is " << gradNorm << ". Successfully quitting." << std::endl);
}

double
LBFGS::norm(const std::vector<double> &x) const
{
    check(x.size() == n_, "");
    double res = 0;
    int i;
#pragma omp parallel for default(shared) private(i) reduction(+:res)
    for(i = 0; i < n_; ++i)
        res += x[i] * x[i];

    return std::sqrt(res);
}

} // namespace Math
