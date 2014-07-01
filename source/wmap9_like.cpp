#include <string>
#include <sstream>
#include <cstring>
//#include <fstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <math_constants.hpp>
#include <wmap9_like.hpp>

#ifdef WMAP9_GFORT
extern "C"
{
    void __wmap_likelihood_9yr_MOD_wmap_likelihood_init();
    void __wmap_likelihood_9yr_MOD_wmap_likelihood_compute(double* cltt, double* clte, double* clee, double* clbb, double* like);
}

extern "C" int __wmap_options_MOD_ttmax;
extern "C" int __wmap_options_MOD_ttmin;
extern "C" int __wmap_options_MOD_lowl_max;
extern "C" bool __wmap_options_MOD_use_lowl_tt;
extern "C" bool __wmap_options_MOD_use_lowl_pol;
extern "C" bool __wmap_options_MOD_use_tt;
extern "C" bool __wmap_options_MOD_use_te;
extern "C" bool __wmap_options_MOD_use_gibbs;
extern "C" bool __wmap_options_MOD_use_tt_beam_ptsrc;
#else
extern "C"
{
    void wmap_likelihood_9yr_mp_wmap_likelihood_init_();
    void wmap_likelihood_9yr_mp_wmap_likelihood_compute_(double* cltt, double* clte, double* clee, double* clbb, double* like);
}

extern "C" int wmap_options_mp_ttmax_;
extern "C" int wmap_options_mp_ttmin_;
extern "C" int wmap_options_mp_lowl_max_;
extern "C" bool wmap_options_mp_use_lowl_tt_;
extern "C" bool wmap_options_mp_use_lowl_pol_;
extern "C" bool wmap_options_mp_use_tt_;
extern "C" bool wmap_options_mp_use_te_;
extern "C" bool wmap_options_mp_use_gibbs_;
extern "C" bool wmap_options_mp_use_tt_beam_ptsrc_;
#endif

bool WMAP9Likelihood::initialized_ = false;

WMAP9Likelihood::WMAP9Likelihood(bool useLowlT, bool useHighlT, bool useLowlP, bool useHighlP, bool useGibbs, bool useTTBeam) : useLowlT_(useLowlT), useHighlT_(useHighlT), useLowlP_(useLowlP), useHighlP_(useHighlP), useGibbs_(useGibbs), useTTBeam_(useTTBeam), like_(10, 0.0), cosmoParams_(6), prevCosmoCalculated_(false)
{
    check(!initialized_, "WMAP 9 likelihood can be initialized only once through the runtime of the program");

#ifdef WMAP9_GFORT
    const int lMin = __wmap_options_MOD_ttmin;
    int lMax = __wmap_options_MOD_ttmax;
#else
    const int lMin = wmap_options_mp_ttmin_;
    int lMax = wmap_options_mp_ttmax_;
#endif
    output_screen1("lMax = " << lMax << std::endl);
    lMax += 1000;

#ifdef WMAP9_GFORT
    __wmap_options_MOD_use_lowl_tt = useLowlT;
    __wmap_options_MOD_use_tt = useHighlT;
    __wmap_options_MOD_use_lowl_pol = useLowlP;
    __wmap_options_MOD_use_te = useHighlP;
    __wmap_options_MOD_use_gibbs = useGibbs;
    __wmap_options_MOD_lowl_max = (useGibbs ? 32 : 30);
    __wmap_options_MOD_use_tt_beam_ptsrc = useTTBeam;
#else
    wmap_options_mp_use_lowl_tt_ = useLowlT;
    wmap_options_mp_use_tt_ = useHighlT;
    wmap_options_mp_use_lowl_pol_ = useLowlP;
    wmap_options_mp_use_te_ = useHighlP;
    wmap_options_mp_use_gibbs_ = useGibbs;
    wmap_options_mp_lowl_max_ = (useGibbs ? 32 : 30);
    wmap_options_mp_use_tt_beam_ptsrc_ = useTTBeam;
#endif

    cmb_.preInitialize(lMax);

#ifdef WMAP9_GFORT
    __wmap_likelihood_9yr_MOD_wmap_likelihood_init();
#else
    wmap_likelihood_9yr_mp_wmap_likelihood_init_();
#endif

    clTT_.resize(lMax + 1, 0.0);
    clTE_.resize(lMax + 1, 0.0);
    clEE_.resize(lMax + 1, 0.0);
    clTB_.resize(lMax + 1, 0.0);
    clEB_.resize(lMax + 1, 0.0);
    clBB_.resize(lMax + 1, 0.0);

    highlT_ = 0;
    lowlTChi2_ = 0;
    lowlTDet_ = 0;
    ttBeam_ = 0;
    teChi2_ = 0;
    teDet_ = 0;
    lowlPChi2_ = 0;
    lowlPDet_ = 0;
    tbChi2_ = 0;
    tbDet_ = 0;

    initialized_ = true;
}

void
WMAP9Likelihood::setCosmoParams(const CosmologicalParams& params)
{
    const bool wantT = true;
    const bool wantPol = (useLowlP_ || useHighlP_);
    const bool wantLens = true;

    cmb_.initialize(params, wantT, wantPol, wantLens);
}

void
WMAP9Likelihood::calculateCls()
{
    const bool wantT = true;
    const bool wantPol = (useLowlP_ || useHighlP_);

    std::vector<double>* tt = &clTT_;
    std::vector<double>* ee = (wantPol ? &clEE_ : NULL);
    std::vector<double>* te = (wantPol ? &clTE_ : NULL);
    std::vector<double>* bb = (wantPol ? &clBB_ : NULL);
    
    cmb_.getLensedCl(tt, ee, te, bb);

    for(int l = 0; l < clTT_.size(); ++l)
        clTT_[l] *= (l * (l + 1) / (2 * Math::pi));

    for(int l = 0; l < clTE_.size(); ++l)
        clTE_[l] *= (l * (l + 1) / (2 * Math::pi));

    for(int l = 0; l < clEE_.size(); ++l)
        clEE_[l] *= (l * (l + 1) / (2 * Math::pi));

    for(int l = 0; l < clBB_.size(); ++l)
        clBB_[l] *= (l * (l + 1) / (2 * Math::pi));

    /*
    std::ofstream out("test_cls.dat");
    for(int l = 2; l < clTT_.size(); ++l)
        out << l << " " << clTT_[l] << " " << clEE_[l] << " " << clBB_[l] << " " << clTE_[l] << " " << 0.0 << std::endl;

    out.close();
    */
}

double
WMAP9Likelihood::calculate(double* params, int nPar)
{
    check(nPar == 6, "");

    const double pivot = 0.05;

    bool needToCalculate = true;
    if(prevCosmoCalculated_)
    {
        check(cosmoParams_.size() == 6, "");
        bool areParamsEqual = true;
        for(int i = 0; i < 6; ++i)
        {
            if(params[i] != cosmoParams_[i])
            {
                areParamsEqual = false;
                break;
            }
        }

        needToCalculate = !areParamsEqual;
    }

    if(needToCalculate)
    {
        LambdaCDMParams lcdmParams(params[0], params[1], params[2], params[3], params[4], std::exp(params[5]) / 1e10, pivot);
        setCosmoParams(lcdmParams);
        calculateCls();

        for(int i = 0; i < 6; ++i)
            cosmoParams_[i] = params[i];

        prevCosmoCalculated_ = true;
    }

    return likelihood();
}

double
WMAP9Likelihood::likelihood()
{
    for(int i = 0; i < like_.size(); ++i)
        like_[i] = 0;

#ifdef WMAP9_GFORT
    __wmap_likelihood_9yr_MOD_wmap_likelihood_compute(&(clTT_[2]), &(clTE_[2]), &(clEE_[2]), &(clBB_[2]), &(like_[0]));
#else
    wmap_likelihood_9yr_mp_wmap_likelihood_compute_(&(clTT_[2]), &(clTE_[2]), &(clEE_[2]), &(clBB_[2]), &(like_[0]));
#endif

    double l = 0;
    for(int i = 0; i < like_.size(); ++i)
        l += 2 * like_[i];

    highlT_ = 2 * like_[0];
    lowlTChi2_ = 2 * like_[1];
    lowlTDet_ = 2 * like_[2];
    ttBeam_ = 2 * like_[3];
    teChi2_ = 2 * like_[4];
    teDet_ = 2 * like_[5];
    lowlPChi2_ = 2 * like_[6];
    lowlPDet_ = 2 * like_[7];
    tbChi2_ = 2 * like_[8];
    tbDet_ = 2 * like_[9];

    return l;
}

WMAP9Likelihood::~WMAP9Likelihood()
{
}

