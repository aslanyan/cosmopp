#ifndef COSMO_PP_LIKELIHOOD_FUNCTION_HPP
#define COSMO_PP_LIKELIHOOD_FUNCTION_HPP

#include <cosmological_params.hpp>

namespace Math
{

/// An abstract likelihood function class
class LikelihoodFunction
{
public:
    virtual ~LikelihoodFunction() {}

    /// A function that calculates the likelihood from the given parameters (purely virtual).
    /// \param params The parameters vector (passed as a pointer to the first element).
    /// \param nParams The number of the parameters.
    /// \return -2ln(likelihood).
    virtual double calculate(double* params, int nParams) = 0;

    /// This function can be used instead of calculate if calculate uses some approximation. In that case calculateExact
    /// can be implemented to return the exact result. The default implementation just calles calculate.
    /// \param params The parameters vector (passed as a pointer to the first element).
    /// \param nParams The number of the parameters.
    /// \return -2ln(likelihood).
    virtual double calculateExact(double* params, int nParams)
    {
        return calculate(params, nParams);
    }
};

class LikelihoodWithDerivs : public LikelihoodFunction
{
public:
    virtual ~LikelihoodWithDerivs() {}

    /// A function that calculates the likelihood from the given parameters (purely virtual).
    /// \param params The parameters vector (passed as a pointer to the first element).
    /// \param nParams The number of the parameters.
    /// \return -2ln(likelihood).
    virtual double calculate(double* params, int nParams) = 0;

    /// A function that calculates the partial derivative of the likelihood from the given parameters (purely virtual).
    /// \param params The parameters vector (passed as a pointer to the first element).
    /// \param nParams The number of the parameters.
    /// \param i The index of the parameter with respect to which the derivative is taken. Should be between 0 and nParams - 1.
    /// \return -2 partial deriv of ln(likelihood) wrt param i.
    virtual double calculateDeriv(double *params, int nParams, int i) = 0;
};

class CosmoLikelihood : public LikelihoodFunction
{
public:
    virtual ~CosmoLikelihood() {}

    /// A function that calculates the likelihood from the given parameters (purely virtual).
    /// \param params The parameters vector (passed as a pointer to the first element).
    /// \param nParams The number of the parameters.
    /// \return -2ln(likelihood).
    virtual double calculate(double* params, int nParams) = 0;

    /// This function sets the cosmological parameters member object of the
    /// likelihood object to be equal to params and also initializes the Class
    /// object. Should be called before calling likelihood().
    /// Example call: setCosmoParams(params)
    /// \param params The cosmological parameters object.
    virtual void setCosmoParams(const CosmologicalParams& params) = 0;

    /// This function calculates the likelihood after setting the cosmological
    /// parameters through setCosmoParams().
    virtual double likelihood() = 0;

    /// This function sets the cosmological model parameters member object of
    /// the likelihood function to be equal to params. It also sets the member
    /// vector of cosmological parameters for use with the calculate function.
    /// Should be called before running one of the parameter sampling
    /// algorithm, i.e. before calling calculate.
    /// Example call: setModelCosmoParams(&params)
    /// \params A pointer to the model cosmological parameters object
    virtual void setModelCosmoParams(CosmologicalParams *params) = 0;
};

} // namespace Math

#endif

