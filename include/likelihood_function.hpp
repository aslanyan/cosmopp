#ifndef COSMO_CPP_LIKELIHOOD_FUNCTION_HPP
#define COSMO_CPP_LIKELIHOOD_FUNCTION_HPP

namespace Math
{

class LikelihoodFunction
{
public:
    virtual ~LikelihoodFunction() {}

    virtual double calculate(double* params, int nParams) = 0;
};

} // namespace Math

#endif

