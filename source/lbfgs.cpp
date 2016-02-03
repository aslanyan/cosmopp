#include <cmath>
#include <algorithm>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <lbfgs.hpp>

/*
#ifdef COSMO_OMP
#include <omp.h>
#endif
*/

namespace Math
{

LBFGS::LBFGS(int n, const RealFunctionMultiDim& f, const RealFunctionMultiToMulti& grad, int m) : n_(n), f_(f), grad_(grad), m_(m), mpi_(CosmoMPI::create())
{
    check(n_ > 0, "the function should have a positive number of parameters");
    check(m_ > 0, "");

    alphaTag_ = mpi_.getCommTag();
    betaTag_ = mpi_.getCommTag();
    z2gTag_ = mpi_.getCommTag();
    stopTag_ = mpi_.getCommTag();
    convergedTag_ = mpi_.getCommTag();
    gradConvergedTag_ = mpi_.getCommTag();
    ysTag_ = mpi_.getCommTag();
    yyTag_ = mpi_.getCommTag();
}

LBFGS::~LBFGS()
{
}


double
LBFGS::minimize(std::vector<double> *x, double epsilon, double gNormTol, int maxIter, void (*callback)(int, double, double, const std::vector<double>&))
{
    mpi_.barrier();

    check(x->size() == n_, "");
    check(epsilon > 0, "");
    check(gNormTol >= 0, "");

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
        callback(iter, val, gradNorm, *x);

    while(true)
    {
        mpi_.barrier();

        q = g;

        const int m = std::min(m_, iter); // use this many previous things
        for(int i = 0; i < m; ++i)
        {
            double dotProduct = 0;
            for(int j = 0; j < n_; ++j)
                dotProduct += s[i][j] * q[j];
            double totalDotProduct = 0;
#ifdef COSMO_MPI
            mpi_.reduce(&dotProduct, &totalDotProduct, 1, CosmoMPI::DOUBLE, CosmoMPI::SUM);
#else
            totalDotProduct = dotProduct;
#endif
            if(mpi_.isMaster())
            {
                alpha[i] = rho[i] * totalDotProduct;
                for(int k = 1; k < mpi_.numProcesses(); ++k)
                    mpi_.send(k, &(alpha[i]), 1, CosmoMPI::DOUBLE, alphaTag_ + k);
            }
            else
            {
                const int k = mpi_.processId();
                check(k != 0, "");
                mpi_.recv(0, &(alpha[i]), 1, CosmoMPI::DOUBLE, alphaTag_ + k);
            }

            for(int j = 0; j < n_; ++j)
                q[j] -= alpha[i] * y[i][j];
        }

        for(int i = 0; i < n_; ++i)
            z[i] = H0k * q[i];

        for(int i = m - 1; i >= 0; --i)
        {
            double dotProduct = 0;
            for(int j = 0; j < n_; ++j)
                dotProduct += y[i][j] * z[j];
            double totalDotProduct = 0;
#ifdef COSMO_MPI
            mpi_.reduce(&dotProduct, &totalDotProduct, 1, CosmoMPI::DOUBLE, CosmoMPI::SUM);
#else
            totalDotProduct = dotProduct;
#endif
            double beta;
            if(mpi_.isMaster())
            {
                beta = rho[i] * totalDotProduct;
                for(int k = 1; k < mpi_.numProcesses(); ++k)
                    mpi_.send(k, &beta, 1, CosmoMPI::DOUBLE, betaTag_ + k);
            }
            else
            {
                const int k = mpi_.processId();
                check(k != 0, "");
                mpi_.recv(0, &beta, 1, CosmoMPI::DOUBLE, betaTag_ + k);
            }

            for(int j = 0; j < n_; ++j)
                z[j] += (alpha[i] - beta) * s[i][j];
        }

        double rate = 1.0;
        double zg = 0;
        for(int i = 0; i < n_; ++i)
            zg += z[i] * g[i];

        double totalZG = 0;
#ifdef COSMO_MPI
            mpi_.reduce(&zg, &totalZG, 1, CosmoMPI::DOUBLE, CosmoMPI::SUM);
#else
            totalZG = zg;
#endif

        int setZToG = 0;
        if(mpi_.isMaster())
        {
            if(totalZG <= 0)
            {
                output_screen("LBFGS iteration " << iter << ": Weird stuff! The descent direction is not a descent direction. Using conjugate gradient at this step!" << std::endl);
                setZToG = 1;
                totalZG = gradNorm * gradNorm;
            }

            for(int k = 1; k < mpi_.numProcesses(); ++k)
                mpi_.send(k, &setZToG, 1, CosmoMPI::INT, z2gTag_ + k);

        }
        else
        {
            const int k = mpi_.processId();
            check(k != 0, "");
            mpi_.recv(0, &setZToG, 1, CosmoMPI::INT, z2gTag_ + k);
        }

        if(setZToG)
            z = g;

        const double tau = 0.5, c = 0.1;

        for(int i = 0; i < n_; ++i)
            searchX[i] = x->at(i) - rate * z[i];
        double newVal = f_.evaluate(searchX);
        int searchIter = 0;
        while(true)
        {
            mpi_.barrier();
            int stop = 0;
            if(mpi_.isMaster())
            {
                if(val - newVal >= rate * c * totalZG || searchIter > 1000)
                    stop = 1;
                for(int k = 1; k < mpi_.numProcesses(); ++k)
                    mpi_.send(k, &stop, 1, CosmoMPI::INT, stopTag_ + k);
            }
            else
            {
                const int k = mpi_.processId();
                check(k != 0, "");
                mpi_.recv(0, &stop, 1, CosmoMPI::INT, stopTag_ + k);
            }
            if(stop)
                break;

            rate *= tau;

            for(int i = 0; i < n_; ++i)
                searchX[i] = x->at(i) - rate * z[i];
            newVal = f_.evaluate(searchX);
            ++searchIter;
        }

        // now move
        for(int i = 0; i < n_; ++i)
            (*x)[i] -= rate * z[i];
        grad_.evaluate(*x, &g);
        const double oldVal = val;
        val = f_.evaluate(*x);
        check(g.size() == n_, "");
        gradNorm = norm(g);

        const double deltaVal = std::abs(val - oldVal);
        const double valMax = std::max(val, oldVal);
        const double ratio = deltaVal / std::max(valMax, 1.0);
        const int minIter = std::max(10, n_);

        int converged = 0;
        if(mpi_.isMaster())
        {
            if(ratio < epsilon && iter >= minIter)
                converged = 1;
            for(int k = 1; k < mpi_.numProcesses(); ++k)
                mpi_.send(k, &converged, 1, CosmoMPI::INT, convergedTag_ + k);
        }
        else
        {
            const int k = mpi_.processId();
            check(k != 0, "");
            mpi_.recv(0, &converged, 1, CosmoMPI::INT, convergedTag_ + k);
        }

        if(converged)
        {
            if(mpi_.isMaster())
            {
                output_screen("LBFGS has reached the required precision!" << std::endl);
            }
            break;
        }

        int gradConverged = 0;
        if(mpi_.isMaster())
        {
            if(gradNorm < gNormTol)
                gradConverged = 1;
            for(int k = 1; k < mpi_.numProcesses(); ++k)
                mpi_.send(k, &gradConverged, 1, CosmoMPI::INT, gradConvergedTag_ + k);
        }
        else
        {
            const int k = mpi_.processId();
            check(k != 0, "");
            mpi_.recv(0, &gradConverged, 1, CosmoMPI::INT, gradConvergedTag_ + k);
        }

        if(gradConverged)
        {
            if(mpi_.isMaster())
            {
                output_screen("LBFGS gradient norm is now below the required tolerance!" << std::endl);
            }
            break;
        }

        // move everything down
        for(int i = m_ - 1; i > 0; --i)
        {
            std::swap(s[i], s[i - 1]);
            std::swap(y[i], y[i - 1]);
            rho[i] = rho[i - 1];
        }

        // set the 0 element
        double ys = 0, yy = 0;
        for(int i = 0; i < n_; ++i)
        {
            s[0][i] = x->at(i) - xPrev[i];
            y[0][i] = g[i] - gPrev[i];
            ys += s[0][i] * y[0][i];
            yy += y[0][i] * y[0][i];
        }
        double totalYS = 0, totalYY = 0;
#ifdef COSMO_MPI
        mpi_.reduce(&ys, &totalYS, 1, CosmoMPI::DOUBLE, CosmoMPI::SUM);
        mpi_.reduce(&yy, &totalYY, 1, CosmoMPI::DOUBLE, CosmoMPI::SUM);
#else
        totalYS = ys;
        totalYY = yy;
#endif

        if(mpi_.isMaster())
        {
            for(int k = 1; k < mpi_.numProcesses(); ++k)
            {
                mpi_.send(k, &totalYS, 1, CosmoMPI::DOUBLE, ysTag_ + k);
                mpi_.send(k, &totalYY, 1, CosmoMPI::DOUBLE, yyTag_ + k);
            }
        }
        else
        {
            const int k = mpi_.processId();
            check(k != 0, "");
            mpi_.recv(0, &totalYS, 1, CosmoMPI::DOUBLE, ysTag_ + k);
            mpi_.recv(0, &totalYY, 1, CosmoMPI::DOUBLE, yyTag_ + k);
        }

        if(totalYS == 0 || totalYY == 0)
        {
            if(mpi_.isMaster())
            {
                output_screen("LBFGS has not moved for some weird reason! Quitting." << std::endl);
            }

            break;
        }

        rho[0] = 1 / totalYS;

        // set H0k
        H0k = totalYS / totalYY;

        for(int i = 0; i < n_; ++i)
        {
            xPrev[i] = x->at(i);
            gPrev[i] = g[i];
        }
        ++iter;

        if(callback)
            callback(iter, val, gradNorm, *x);

        if(iter > maxIter)
        {
            if(mpi_.isMaster())
            {
                output_screen("LBFGS has reached the maximum number of iterations of " << maxIter << ". Quitting." << std::endl);
            }

            break;
        }
    }

    if(mpi_.isMaster())
    {
        output_screen("LBFGS has converged after " << iter << " iterations. Successfully quitting." << std::endl);
        output_screen("Iterations: " << iter << ", function value: " << val << ", gradient norm: " << gradNorm << std::endl);
    }
    return val;
}

double
LBFGS::norm(const std::vector<double> &x) const
{
    check(x.size() == n_, "");
    double res = 0;
    int i;
//#pragma omp parallel for default(shared) private(i) reduction(+:res)
    for(i = 0; i < n_; ++i)
        res += x[i] * x[i];

    double totalRes = 0;
#ifdef COSMO_MPI
    mpi_.reduce(&res, &totalRes, 1, CosmoMPI::DOUBLE, CosmoMPI::SUM);
#else
    totalRes = res;
#endif

    return std::sqrt(totalRes);
}

} // namespace Math
