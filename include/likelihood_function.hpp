#ifndef COSMO_PP_LIKELIHOOD_FUNCTION_HPP
#define COSMO_PP_LIKELIHOOD_FUNCTION_HPP

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

} // namespace Math

#endif

