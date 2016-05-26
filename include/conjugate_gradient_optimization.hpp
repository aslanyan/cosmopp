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

class CG
{
public:
    enum Method {FLETCHER_REEVES = 0, POLAK_RIBIERE, HESTENES_STIEFEL, DAI_YUAN, METHOD_MAX };

    CG(int n, const RealFunctionMultiDim& f, const RealFunctionMultiToMulti& grad, const std::vector<double>& starting) : factory_(n), f_(f, grad)
    {
        s_ = factory_.giveMeOne();
        s_->contents() = starting;
        cg_.reset(new CG_General<BasicLargeVector, BasicLargeVectorFactory, BasicLBFGSFunc>(&factory_, &f_, *s_)); 
    }

    void setStarting(const std::vector<double>& starting)
    {
        s_->contents() = starting;
        cg_->setStarting(*s_);
    }

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

