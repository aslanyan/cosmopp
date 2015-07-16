#include <macros.hpp>
#include <modecode.hpp>

#ifdef MODECODE_GFORT
extern "C"
{
    void __access_modpk_MOD_potinit();
    void __access_modpk_MOD_evolve(double*, double*, double*);
}

extern "C" bool __camb_interface_MOD_modpkoutput;
extern "C" bool __modpkparams_MOD_modpk_physical_priors;
extern "C" int __modpkparams_MOD_potential_choice;
extern "C" bool __modpkparams_MOD_flag_do_reconstruction;
extern "C" bool __modpkparams_MOD_vnderivs;
extern "C" bool __modpkparams_MOD_instreheat;
extern "C" double __modpkparams_MOD_k_pivot;
extern "C" double __modpkparams_MOD_n_pivot;
extern "C" bool __modpkparams_MOD_slowroll_infl_end;
extern "C" double __modpkparams_MOD_vparams;
extern "C" double __modpkparams_MOD_findiffdphi;
extern "C" int __camb_interface_MOD_pk_bad;
extern "C" bool __modpkparams_MOD_eternal_infl_ok;
extern "C" double __modpkparams_MOD_modpk_w_primordial_lower;
extern "C" double __modpkparams_MOD_modpk_w_primordial_upper;
extern "C" double __modpkparams_MOD_modpk_rho_reheat;
#else
extern "C"
{
    void access_modpk_mp_potinit_();
    void access_modpk_mp_evolve_(double*, double*, double*);
}

extern "C" bool camb_interface_mp_modpkoutput_;
extern "C" bool modpkparams_mp_modpk_physical_priors_;
extern "C" int modpkparams_mp_potential_choice_;
extern "C" bool modpkparams_mp_flag_do_reconstruction_;
extern "C" bool modpkparams_mp_vnderivs_;
extern "C" bool modpkparams_mp_instreheat_;
extern "C" double modpkparams_mp_k_pivot_;
extern "C" double modpkparams_mp_n_pivot_;
extern "C" bool modpkparams_mp_slowroll_infl_end_;
extern "C" double modpkparams_mp_vparams_;
extern "C" double modpkparams_mp_findiffdphi_;
extern "C" int camb_interface_mp_pk_bad_;
extern "C" bool modpkparams_mp_eternal_infl_ok_;
extern "C" double modpkparams_mp_modpk_w_primordial_lower_;
extern "C" double modpkparams_mp_modpk_w_primordial_upper_;
extern "C" double modpkparams_mp_modpk_rho_reheat_;
#endif

int ModeCode::nVPar_ = 0;
double* ModeCode::vParams_ = NULL;
Math::TableFunction<double, double> ModeCode::scalarPs_;
Math::TableFunction<double, double> ModeCode::tensorPs_;
std::map<double, double> ModeCode::scalarLower_;
std::map<double, double> ModeCode::scalarUpper_;
std::map<double, double> ModeCode::tensorLower_;
std::map<double, double> ModeCode::tensorUpper_;


void
ModeCode::initialize(int potentialChoice, double kPivot, double NPivot, bool instantReheating, bool physicalPriors, bool slowRollEnd, bool eternalInflOK, double kMin, double kMax, int nPoints)
{
    check(potentialChoice >= 1, "invalid potential choice " << potentialChoice);
    check(kPivot > 0, "invalid k_pivot = " << kPivot);
    check(NPivot > 0, "invalid N_pivot = " << NPivot);
    check(kMin > 0, "invalid k_min = " << kMin);
    check(kMax >= kMin, "invalid k_max = " << kMax);
    check(nPoints > 0, "invalid nPoints = " << nPoints << ", needs to be positive");

    bool outputFlag = false;
#ifdef VERBOSE2
    outputFlag = true;
#endif

#ifdef MODECODE_GFORT
    __camb_interface_MOD_modpkoutput = outputFlag;
    __modpkparams_MOD_modpk_physical_priors = physicalPriors;
    __modpkparams_MOD_modpk_rho_reheat = 1e44;
    __modpkparams_MOD_modpk_w_primordial_lower = -1.0 / 3.0;
    __modpkparams_MOD_modpk_w_primordial_upper = 1.0;
    __modpkparams_MOD_potential_choice = potentialChoice;
    __modpkparams_MOD_flag_do_reconstruction = false;
    __modpkparams_MOD_vnderivs = false;
    __modpkparams_MOD_instreheat = instantReheating;
    __modpkparams_MOD_k_pivot = kPivot;
    __modpkparams_MOD_n_pivot = NPivot;
    __modpkparams_MOD_slowroll_infl_end = slowRollEnd;
    __modpkparams_MOD_findiffdphi = 1.0e-16;
    __modpkparams_MOD_eternal_infl_ok = eternalInflOK;

    vParams_ = &__modpkparams_MOD_vparams;
#else
    camb_interface_mp_modpkoutput_ = outputFlag;
    modpkparams_mp_modpk_physical_priors_ = physicalPriors;
    modpkparams_mp_modpk_rho_reheat_ = 1e44;
    modpkparams_mp_modpk_w_primordial_lower_ = -1.0 / 3.0;
    modpkparams_mp_modpk_w_primordial_upper_ = 1.0;
    modpkparams_mp_potential_choice_ = potentialChoice;
    modpkparams_mp_flag_do_reconstruction_ = false;
    modpkparams_mp_vnderivs_ = false;
    modpkparams_mp_instreheat_ = instantReheating;
    modpkparams_mp_k_pivot_ = kPivot;
    modpkparams_mp_n_pivot_ = NPivot;
    modpkparams_mp_slowroll_infl_end_ = slowRollEnd;
    modpkparams_mp_findiffdphi_ = 1.0e-16;
    modpkparams_mp_eternal_infl_ok_ = eternalInflOK;

    vParams_ = &modpkparams_mp_vparams_;
#endif

    const double deltaK = (kMax - kMin) / nPoints;

    scalarPs_.clear();
    tensorPs_.clear();

    for(int i = 0; i <= nPoints; ++i)
    {
        const double k = (i == nPoints ? kMax : kMin + i * deltaK);
        scalarPs_[k] = 0;
        tensorPs_[k] = 0;
        scalarLower_[k] = 1e-10;
        scalarUpper_[k] = 1e-8;
        tensorLower_[k] = 0;
        tensorUpper_[k] = 1e-8;
    }

    switch(potentialChoice)
    {
    case 1:
        nVPar_ = 1;
        break;
    case 12:
        nVPar_ = 5;
        break;
    default:
        check(false, "not implemented, TBD better");
        break;
    }
}

void
ModeCode::setNPivot(double NPivot)
{
    check(NPivot > 0, "invalid N_pivot = " << NPivot);
#ifdef MODECODE_GFORT
    __modpkparams_MOD_n_pivot = NPivot;
#else
    modpkparams_mp_n_pivot_ = NPivot;
#endif
}

double
ModeCode::getNPivot()
{
    double nPiv;
#ifdef MODECODE_GFORT
    nPiv = __modpkparams_MOD_n_pivot;
#else
    nPiv = modpkparams_mp_n_pivot_;
#endif
    return nPiv;
}

void
ModeCode::addKValue(double k, double sMin, double sMax, double tMin, double tMax)
{
    check(k > 0, "");
    check(scalarPs_.find(k) == scalarPs_.end(), "already exists");
    scalarPs_[k] = 0;
    tensorPs_[k] = 0;
    scalarLower_[k] = sMin;
    scalarUpper_[k] = sMax;
    tensorLower_[k] = tMin;
    tensorUpper_[k] = tMax;
}

bool
ModeCode::calculate(const std::vector<double>& vParams)
{
    check(vParams.size() <= 10, "too many params");

    for(int i = 0; i < vParams.size(); ++i)
        vParams_[i] = vParams[i];

#ifdef MODECODE_GFORT
    __access_modpk_MOD_potinit();
#else
    access_modpk_mp_potinit_();
#endif

#ifdef MODECODE_GFORT
    if(__camb_interface_MOD_pk_bad != 0)
        return false;
#else
    if(camb_interface_mp_pk_bad_ != 0)
        return false;
#endif

    Math::TableFunction<double, double>::iterator it = scalarPs_.begin(), it1 = tensorPs_.begin();
    for(; it != scalarPs_.end(); ++it)
    {
        check(it1 != tensorPs_.end(), "");
        double k = (*it).first;
        check((*it1).first == k, "");
        double s, t;
        check(scalarLower_.find(k) != scalarLower_.end(), "");
        check(scalarUpper_.find(k) != scalarUpper_.end(), "");
        check(tensorLower_.find(k) != tensorLower_.end(), "");
        check(tensorUpper_.find(k) != tensorUpper_.end(), "");
        const double sMin = scalarLower_[k], sMax = scalarUpper_[k], tMin = tensorLower_[k], tMax = tensorUpper_[k];
#ifdef MODECODE_GFORT
        if(__camb_interface_MOD_pk_bad != 0)
            return false;
        __access_modpk_MOD_evolve(&k, &s, &t);
        if(__camb_interface_MOD_pk_bad != 0 || s <= sMin || s >= sMax || t <= tMin || t >= tMax)
            return false;
#else
        if(camb_interface_mp_pk_bad_ != 0)
            return false;
        access_modpk_mp_evolve_(&k, &s, &t);
        if(camb_interface_mp_pk_bad_ != 0 || s <= sMin || s >= sMax || t <= tMin || t >= tMax)
            return false;
#endif
        (*it).second = s;
        (*it1).second = t;
        ++it1;
    }

    return true;
}
