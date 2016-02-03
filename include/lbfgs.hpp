#ifndef COSMO_PP_LBFGS_HPP
#define COSMO_PP_LBFGS_HPP

#include <function.hpp>
#include <cosmo_mpi.hpp>

namespace Math
{

class LBFGS
{
public:
    LBFGS(int n, const RealFunctionMultiDim& f, const RealFunctionMultiToMulti& grad, int m = 10);
    ~LBFGS();

    double minimize(std::vector<double> *x, double epsilon = 1e-3, double gNormTol = 1e-5, int maxIter = 1000000, void (*callback)(int, double, double, const std::vector<double>&) = NULL);

private:
    double norm(const std::vector<double>& x) const;

private:
    const RealFunctionMultiDim& f_;
    const RealFunctionMultiToMulti& grad_;
    const int n_;
    const int m_;

    CosmoMPI& mpi_;

    int alphaTag_, betaTag_, z2gTag_, stopTag_, convergedTag_, gradConvergedTag_, ysTag_, yyTag_;
};

} // namespace Math

#endif

