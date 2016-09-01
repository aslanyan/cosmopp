#ifndef COSMO_PP_CG_OPTIMIZATION_HPP
#define COSMO_PP_CG_OPTIMIZATION_HPP

#include <memory>

#include <macros.hpp>
#include <function.hpp>
#include <cosmo_mpi.hpp>
#include <conjugate_gradient_general.hpp>
#include <lbfgs.hpp>

namespace Math
{


/// A Non-linear conjugate gradient optimizer. This is a simpler version which would be used for many applications. For a more general interface use CG_GENERAL.
class CG
{
public:
    /// Choice of the nonlinear method.
    enum Method {FLETCHER_REEVES = 0, POLAK_RIBIERE, HESTENES_STIEFEL, DAI_YUAN, METHOD_MAX };

    /// Constructor.
    /// \param n The number of dimensions.
    /// \param f The function to minimize.
    /// \param grad The gradient of the function.
    /// \param starting The starting point.
    CG(int n, const RealFunctionMultiDim& f, const RealFunctionMultiToMulti& grad, const std::vector<double>& starting) : factory_(n), f_(f, grad)
    {
        s_ = factory_.giveMeOne();
        s_->contents() = starting;
        cg_.reset(new CG_General<BasicLargeVector, BasicLargeVectorFactory, BasicLBFGSFunc>(&factory_, &f_, *s_)); 
    }

    /// Set a new starting point. This will completely reset the minimizer so you can run minimize again from the new starting point.
    /// \param starting The new starting point.
    void setStarting(const std::vector<double>& starting)
    {
        s_->contents() = starting;
        cg_->setStarting(*s_);
    }

    /// Function for minimization (the main function of this class).
    /// \param res A pointer to a vector where the result will be stored.
    /// \param epsilon Threshold for optimization. Will stop if two successive iterations change the function value by less than epsilon.
    /// \param gNormTol Threshold for the norm of the gradient. The minimizer will stop if the gradient norm at the given iteration is less than gNormTol.
    /// \param maxIter Maximum number of iterations. The optimizer will stop if it reaches this maximum number, regardless of convergence.
    /// \param m The nonlinear method to use.
    /// \param callback A callback function to be called at each iteration. If NULL (default) then no callback will be used.
    double minimize(std::vector<double> *res, double epsilon = 1e-3, double gNormTol = 1e-5, int maxIter = 1000000, Method m = FLETCHER_REEVES, void (*callback)(int, double, double, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&) = NULL)
    {
        CG_General<BasicLargeVector, BasicLargeVectorFactory, BasicLBFGSFunc>::Method method;
        switch(m)
        {
        case FLETCHER_REEVES:
            method = CG_General<BasicLargeVector, BasicLargeVectorFactory, BasicLBFGSFunc>::FLETCHER_REEVES;
            break;
        case POLAK_RIBIERE:
            method = CG_General<BasicLargeVector, BasicLargeVectorFactory, BasicLBFGSFunc>::POLAK_RIBIERE;
            break;
        case HESTENES_STIEFEL:
            method = CG_General<BasicLargeVector, BasicLargeVectorFactory, BasicLBFGSFunc>::HESTENES_STIEFEL;
            break;
        case DAI_YUAN:
            method = CG_General<BasicLargeVector, BasicLargeVectorFactory, BasicLBFGSFunc>::DAI_YUAN;
            break;
        default:
            check(false, "invalid method");
            break;
        }
        cb_.set(callback);
        const double val = cg_->minimize(s_, epsilon, gNormTol, maxIter, method, &cb_);
        (*res) = s_->contents();
        return val;
    }

private:
    BasicLargeVectorFactory factory_;
    BasicLBFGSFunc f_;

    BasicLargeVector *s_;

    std::unique_ptr<CG_General<BasicLargeVector, BasicLargeVectorFactory, BasicLBFGSFunc> > cg_;
    BasicLBFGSCallBack cb_;
};

} // namespace Math

#endif

