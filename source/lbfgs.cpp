#include <cmath>
#include <algorithm>
#include <vector>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <lbfgs.hpp>
#include <random.hpp>
#include <cosmo_mpi.hpp>

/*
#ifdef COSMO_OMP
#include <omp.h>
#endif
*/

namespace Math
{

/*
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
    int functionEval = 0;

    double val = f_.evaluate(*x);
    std::vector<double> g;
    grad_.evaluate(*x, &g);
    check(g.size() == n_, "");
    double gradNorm = norm(g);
    ++functionEval;

    std::vector<double> xPrev(*x), gPrev(g);
    std::vector<double> q(n_);
    std::vector<double> z(n_);

    std::vector<double> searchX(n_);

    double H0k = 1;
    if(callback)
        callback(iter, val, gradNorm, *x);

    while(true)
    {
        q = g;

        const int m = std::min(m_, iter); // use this many previous things
        for(int i = 0; i < m; ++i)
        {
            double dotProduct = 0;
            for(int j = 0; j < n_; ++j)
                dotProduct += s[i][j] * q[j];
            double totalDotProduct = dotProduct;
#ifdef COSMO_MPI
            mpi_.reduce(&dotProduct, &totalDotProduct, 1, CosmoMPI::DOUBLE, CosmoMPI::SUM);
#endif
            alpha[i] = rho[i] * totalDotProduct;
            mpi_.bcast(&(alpha[i]), 1, CosmoMPI::DOUBLE);

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
            double totalDotProduct = dotProduct;
#ifdef COSMO_MPI
            mpi_.reduce(&dotProduct, &totalDotProduct, 1, CosmoMPI::DOUBLE, CosmoMPI::SUM);
#endif
            double beta = rho[i] * totalDotProduct;
            mpi_.bcast(&beta, 1, CosmoMPI::DOUBLE);

            for(int j = 0; j < n_; ++j)
                z[j] += (alpha[i] - beta) * s[i][j];
        }

        double zg = 0;
        for(int i = 0; i < n_; ++i)
            zg += z[i] * g[i];

        double totalZG = zg;
#ifdef COSMO_MPI
            mpi_.reduce(&zg, &totalZG, 1, CosmoMPI::DOUBLE, CosmoMPI::SUM);
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
        }
        mpi_.bcast(&setZToG, 1, CosmoMPI::INT);

        if(setZToG)
            z = g;

        const double tau = 0.5, c = 0.01;
        double rate = 1.0;
        for(int i = 0; i < n_; ++i)
            searchX[i] = x->at(i) - rate * z[i];
        double newVal = f_.evaluate(searchX);
        ++functionEval;
        int searchIter = 0;
        while(true)
        {
            mpi_.barrier();
            int stop = 0;
            if(mpi_.isMaster())
            {
                if(val - newVal >= rate * c * totalZG || searchIter > 1000)
                    stop = 1;
            }
            mpi_.bcast(&stop, 1, CosmoMPI::INT);
            if(stop)
                break;

            rate *= tau;

            for(int i = 0; i < n_; ++i)
                searchX[i] = x->at(i) - rate * z[i];
            newVal = f_.evaluate(searchX);
            ++functionEval;
            ++searchIter;
        }

        // now move
        //for(int i = 0; i < n_; ++i)
            //(*x)[i] -= rate * z[i];
        (*x) = searchX;
        const double oldVal = val;
        //val = f_.evaluate(*x);
        val = newVal;
        grad_.evaluate(*x, &g);
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
        }
        mpi_.bcast(&converged, 1, CosmoMPI::INT);

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
        }
        mpi_.bcast(&gradConverged, 1, CosmoMPI::INT);

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
        //double ys = 0, yy = 0;
        double ys_yy[2] = {0, 0};
        for(int i = 0; i < n_; ++i)
        {
            s[0][i] = x->at(i) - xPrev[i];
            y[0][i] = g[i] - gPrev[i];
            ys_yy[0] += s[0][i] * y[0][i];
            ys_yy[1] += y[0][i] * y[0][i];
        }
        //double totalYS = 0, totalYY = 0;
        double total_ys_yy[2];
#ifdef COSMO_MPI
        //mpi_.reduce(&ys, &totalYS, 1, CosmoMPI::DOUBLE, CosmoMPI::SUM);
        //mpi_.reduce(&yy, &totalYY, 1, CosmoMPI::DOUBLE, CosmoMPI::SUM);
        mpi_.reduce(ys_yy, total_ys_yy, 2, CosmoMPI::DOUBLE, CosmoMPI::SUM);
#else
        //totalYS = ys;
        //totalYY = yy;
        total_ys_yy[0] = ys_yy[0];
        total_ys_yy[1] = ys_yy[1];
#endif

        //mpi_.bcast(&totalYS, 1, CosmoMPI::DOUBLE);
        //mpi_.bcast(&totalYY, 1, CosmoMPI::DOUBLE);
        mpi_.bcast(&total_ys_yy, 2, CosmoMPI::DOUBLE);

        if(total_ys_yy[0] == 0 || total_ys_yy[1] == 0)
        {
            if(mpi_.isMaster())
            {
                output_screen("LBFGS has not moved for some weird reason! Quitting." << std::endl);
            }

            break;
        }

        rho[0] = 1 / total_ys_yy[0];

        // set H0k
        H0k = total_ys_yy[0] / total_ys_yy[1];

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
        output_screen("Iterations: " << iter << ", function evaluations: " << functionEval << ", function value: " << val << ", gradient norm: " << gradNorm << std::endl);
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

*/

void
BasicLargeVector::copy(const BasicLargeVector& other, double c)
{
    check(v_.size() == other.v_.size(), "");
    for(int i = 0; i < v_.size(); ++i)
        v_[i] = c * other.v_[i];
}

void
BasicLargeVector::setToZero()
{
    for(int i = 0; i < v_.size(); ++i)
        v_[i] = 0;
}

double
BasicLargeVector::norm() const
{
    return std::sqrt(dotProduct(*this));
}

double
BasicLargeVector::dotProduct(const BasicLargeVector& other) const
{
    check(other.v_.size() == v_.size(), "");
    double s = 0;
    for(int i = 0; i < v_.size(); ++i)
        s += v_[i] * other.v_[i];

    double total = s;
#ifdef COSMO_MPI
    CosmoMPI::create().reduce(&s, &total, 1, CosmoMPI::DOUBLE, CosmoMPI::SUM);
    CosmoMPI::create().bcast(&total, 1, CosmoMPI::DOUBLE);
#endif
    return total;
}

void
BasicLargeVector::add(const BasicLargeVector& other, double c)
{
    check(other.v_.size() == v_.size(), "");
    for(int i = 0; i < v_.size(); ++i)
        v_[i] += c * other.v_[i];
}

void
BasicLargeVector::multiply(const BasicLargeVector& other)
{
    check(other.v_.size() == v_.size(), "");
    for(int i = 0; i < v_.size(); ++i)
        v_[i] *= other.v_[i];
}

void
BasicLargeVector::divide(const BasicLargeVector& other)
{
    check(other.v_.size() == v_.size(), "");
    for(int i = 0; i < v_.size(); ++i)
    {
        check(other.v_[i] != 0, "division by 0 at index" << i);
        v_[i] /= other.v_[i];
    }
}

void
BasicLargeVector::pow(double p)
{
    for(int i = 0; i < v_.size(); ++i)
        v_[i] = std::pow(v_[i], p);
}

void
BasicLargeVector::swap(BasicLargeVector& other)
{
    check(other.v_.size() == v_.size(), "");
    v_.swap(other.v_);
}

void
BasicLBFGSFunc::whitenoise(int seed, BasicLargeVector* x, double amplitude)
{
    const int nProc = CosmoMPI::create().numProcesses();
    const int pid = CosmoMPI::create().processId();
    std::vector<int> seeds(nProc);
    Math::UniformRealGenerator g1(seed, 0, 1);
    for(int i = 0; i < nProc; ++i)
        seeds[i] = static_cast<int>(std::ceil(1000000 * g1.generate()));

    Math::GaussianGenerator g(seeds[pid], 0, 1);

    std::vector<double>& v = x->contents();
    for(int i = 0; i < v.size(); ++i)
        v[i] = amplitude * g.generate();
}

} // namespace Math
