#ifndef COSMO_PP_HMC_GENERAL_HPP
#define COSMO_PP_HMC_GENERAL_HPP

namespace Math
{

/*
class HMCTraits
{
public:
    unsigned long nPar() const;
};
*/

template<typename HMCTraits>
class HMCGeneral
{
public:
    HMCGeneral(HMCTraits *traits);
    ~HMCGeneral();
    void run(int iters);

private:
    HMCTraits *traits_;
};

} // namespace Math

#endif

