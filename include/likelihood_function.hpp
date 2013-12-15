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
};

} // namespace Math

#endif

