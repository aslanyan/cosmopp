#include <string>
#include <sstream>
#include <cstring>
#include <fstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <planck_like.hpp>
#include <timer.hpp>

#include <clik.h>

#define MY_STRINGIZE1(P) #P
#define MY_STRINGIZE(P) MY_STRINGIZE1(P)
#define PLANCK_DATA_DIR_STR MY_STRINGIZE(PLANCK_DATA_DIR)

#ifdef COSMO_PLANCK_15

namespace
{

class PlanckLikelihoodContainer
{
public:
    PlanckLikelihoodContainer() : lowT_(NULL), lowTP_(NULL), highT_(NULL), highTP_(NULL), highTLite_(NULL), highTPLite_(NULL), lensT_(NULL), lensTP_(NULL), planckLikeDir_(PLANCK_DATA_DIR_STR)
    {
    }

    void* getLowT()
    {
        if(lowT_)
            return lowT_;

        std::stringstream path;
        path << planckLikeDir_ << "/low_l/commander/commander_rc2_v1.1_l2_29_B.clik";
        char pathCStr[200];
        std::strcpy(pathCStr, path.str().c_str());
        lowT_ = clik_init(pathCStr, NULL);
        return lowT_;
    }

    void* getLowTP()
    {
        if(lowTP_)
            return lowTP_;

        std::stringstream path;
        path << planckLikeDir_ << "/low_l/bflike/lowl_SMW_70_dx11d_2014_10_03_v5c_Ap.clik";
        char pathCStr[200];
        std::strcpy(pathCStr, path.str().c_str());
        lowTP_ = clik_init(pathCStr, NULL);
        return lowTP_;
    }

    void* getHighT()
    {
        if(highT_)
            return highT_;

        std::stringstream path;
        path << planckLikeDir_ << "/hi_l/plik/plik_dx11dr2_HM_v18_TT.clik";
        char pathCStr[200];
        std::strcpy(pathCStr, path.str().c_str());
        highT_ = clik_init(pathCStr, NULL);
        return highT_;
    }

    void* getHighTP()
    {
        if(highTP_)
            return highTP_;

        std::stringstream path;
        path << planckLikeDir_ << "/hi_l/plik/plik_dx11dr2_HM_v18_TTTEEE.clik";
        char pathCStr[200];
        std::strcpy(pathCStr, path.str().c_str());
        highTP_ = clik_init(pathCStr, NULL);
        return highTP_;
    }

    void* getHighTLite()
    {
        if(highTLite_)
            return highTLite_;

        check(!highTPLite_, "high TP already in use, cannot use high T in addition (it's a problem in Planck likelihood code itself)");

        std::stringstream path;
        path << planckLikeDir_ << "/hi_l/plik_lite/plik_lite_v18_TT.clik";
        char pathCStr[200];
        std::strcpy(pathCStr, path.str().c_str());
        highTLite_ = clik_init(pathCStr, NULL);
        return highTLite_;
    }

    void* getHighTPLite()
    {
        if(highTPLite_)
            return highTPLite_;

        check(!highTLite_, "high T already in use, cannot use high TP in addition (it's a problem in Planck likelihood code itself)");

        std::stringstream path;
        path << planckLikeDir_ << "/hi_l/plik_lite/plik_lite_v18_TTTEEE.clik";
        char pathCStr[200];
        std::strcpy(pathCStr, path.str().c_str());
        highTPLite_ = clik_init(pathCStr, NULL);
        return highTPLite_;
    }

    void* getLensT()
    {
        if(lensT_)
            return lensT_;

        std::stringstream path;
        path << planckLikeDir_ << "/lensing/smica_g30_ftl_full_pttptt.clik_lensing";
        char pathCStr[200];
        std::strcpy(pathCStr, path.str().c_str());
        lensT_ = clik_lensing_init(pathCStr, NULL);
        return lensT_;
    }

    void* getLensTP()
    {
        if(lensTP_)
            return lensTP_;

        std::stringstream path;
        path << planckLikeDir_ << "/lensing/smica_g30_ftl_full_pp.clik_lensing";
        char pathCStr[200];
        std::strcpy(pathCStr, path.str().c_str());
        lensTP_ = clik_lensing_init(pathCStr, NULL);
        return lensTP_;
    }

    ~PlanckLikelihoodContainer()
    {
        if(lowT_) clik_cleanup(&lowT_);
        if(lowTP_) clik_cleanup(&lowTP_);
        if(highT_) clik_cleanup(&highT_);
        if(highTP_) clik_cleanup(&highTP_);
        if(highTLite_) clik_cleanup(&highTLite_);
        if(highTPLite_) clik_cleanup(&highTPLite_);
        if(lensT_) clik_lensing_cleanup(&lensT_);
        if(lensTP_) clik_lensing_cleanup(&lensTP_);
    }

private:
    clik_object *lowT_;
    clik_object *lowTP_;
    clik_object *highT_;
    clik_object *highTP_;
    clik_object *highTLite_;
    clik_object *highTPLite_;
    clik_lensing_object *lensT_;
    clik_lensing_object *lensTP_;
    std::string planckLikeDir_;
};

}

PlanckLikelihood::PlanckLikelihood(bool lowT, bool lowP, bool highT, bool highP, bool highLikeLite, bool lensingT, bool lensingP, bool includeTensors, double kPerDecade, bool useOwnCmb) : spectraNames_(6), lensSpectraNames_(7), low_(NULL), high_(NULL), lens_(NULL), lowT_(lowT), lowP_(lowP), highT_(highT), highP_(highP), highLikeLite_(highLikeLite), lensingT_(lensingT), lensingP_(lensingP), cmb_(NULL), modelParams_(NULL), aPlanck_(1), aPol_(1), szPrior_(false)
{
    check(!lowP || lowT, "cannot include lowP without lowT");
    check(!highP || highT, "cannot include highP without highT");
    check(!lensingP || lensingT, "cannot include lensingP wihtout lensingT");

    check(lowT || highT || lensingT, "at least one likelihood must be specified");

    std::string planckLikeDir = PLANCK_DATA_DIR_STR;
    output_screen("Planck data dir = " << planckLikeDir << std::endl);

    spectraNames_[0] = std::string("TT");
    spectraNames_[1] = std::string("EE");
    spectraNames_[2] = std::string("BB");
    spectraNames_[3] = std::string("TE");
    spectraNames_[4] = std::string("TB");
    spectraNames_[5] = std::string("EB");

    lensSpectraNames_[0] = std::string("PP");
    lensSpectraNames_[1] = std::string("TT");
    lensSpectraNames_[2] = std::string("EE");
    lensSpectraNames_[3] = std::string("BB");
    lensSpectraNames_[4] = std::string("TE");
    lensSpectraNames_[5] = std::string("TB");
    lensSpectraNames_[6] = std::string("EB");

    int lMax[6];
    int lensLMax[7];

    lMax_ = 0;
    lowLMax_ = 0;
    highLMax_ = 0;
    lensLMax_ = 0;

    parname* names;

    static PlanckLikelihoodContainer planckLikelihoodContainer;

    if(lowT)
    {
        if(lowP)
            low_ = planckLikelihoodContainer.getLowTP();
        else
            low_ = planckLikelihoodContainer.getLowT();

        clik_get_lmax(low_, lMax, NULL);
        for(int i = 0; i < 6; ++i)
        {
            output_screen1("Low likelihood has " << spectraNames_[i] << " with l_max = " << lMax[i] << std::endl);
            if(lMax[i] > lowLMax_)
                lowLMax_ = lMax[i];
        }

        if(lowLMax_ > lMax_)
            lMax_ = lowLMax_;

        /*
        int extraParams = clik_get_extra_parameter_names(low_, &names, NULL);
        output_screen1("Low likelihood has " << extraParams << " extra parameters. Their names are as follows:" << std::endl);
        for(int i = 0; i < extraParams; ++i)
        {
            output_screen1(names[i] << std::endl);
        }
        output_screen1(std::endl);
        delete names;
        */
    }

    if(highT)
    {
        if(highLikeLite)
        {
            if(highP)
                high_ = planckLikelihoodContainer.getHighTPLite();
            else
                high_ = planckLikelihoodContainer.getHighTLite();
        }
        else
        {
            if(highP)
                high_ = planckLikelihoodContainer.getHighTP();
            else
                high_ = planckLikelihoodContainer.getHighT();
        }
        clik_get_lmax(high_, lMax, NULL);
        for(int i = 0; i < 6; ++i)
        {
            output_screen1("High likelihood has " << spectraNames_[i] << " with l_max = " << lMax[i] << std::endl);
            if(lMax[i] > highLMax_)
                highLMax_ = lMax[i];
        }

        if(highLMax_ > lMax_)
            lMax_ = highLMax_;

        /*
        int extraParams = clik_get_extra_parameter_names(high_, &names, NULL);
        output_screen1("High likelihood has " << extraParams << " extra parameters. Their names are as follows:" << std::endl);
        for(int i = 0; i < extraParams; ++i)
        {
            output_screen1(names[i] << std::endl);
        }
        output_screen1(std::endl);
        delete names;
        */
    }

    if(lensingT)
    {
        if(lensingP)
            lens_ = planckLikelihoodContainer.getLensTP();
        else
            lens_ = planckLikelihoodContainer.getLensT();

        clik_lensing_get_lmaxs(static_cast<clik_lensing_object*>(lens_), lensLMax, NULL);
        for(int i = 0; i < 7; ++i)
        {
            output_screen1("Lensing likelihood has " << lensSpectraNames_[i] << " with l_max = " << lensLMax[i] << std::endl);
            if(lensLMax[i] > lensLMax_)
                lensLMax_ = lensLMax[i];
        }

        if(lensLMax_ > lMax_)
            lMax_ = lensLMax_;

        /*
        int extraParams = clik_lensing_get_extra_parameter_names(static_cast<clik_lensing_object*>(lens_), &names, NULL);
        output_screen1("Lensing likelihood has " << extraParams << " extra parameters. Their names are as follows:" << std::endl);
        for(int i = 0; i < extraParams; ++i)
        {
            output_screen1(names[i] << std::endl);
        }
        output_screen1(std::endl);
        delete names;
        */
    }

    output_screen1("Low l_max = " << lowLMax_ << std::endl);
    output_screen1("High l_max = " << highLMax_ << std::endl);
    output_screen1("Lensing l_max = " << lensLMax_ << std::endl);
    output_screen1("Total l_max = " << lMax_ << std::endl);

    if(useOwnCmb)
    {
        cmb_ = new CMB;
        cmb_->preInitialize(lMax_ + 1000, false, true, includeTensors, lMax_ + 1000, kPerDecade);
    }

    input_.resize(3 * lMax_ + 10000);
    if(highT)
    {
        highExtra_.resize(highP ? 32 : 15, 0);
        if(highP)
            beamExtra_.resize(60, 0);
    }
}

PlanckLikelihood::~PlanckLikelihood()
{
    if(cmb_)
        delete cmb_;
}

void
PlanckLikelihood::setCosmoParams(const CosmologicalParams& params)
{
    check(cmb_, "own CMB must be used (set in constructor)");

    bool needToCalculate = true;
    params.getAllParameters(currentCosmoParams_);
    if(params.name() == prevCosmoParamsName_)
    {
        check(currentCosmoParams_.size() == prevCosmoParams_.size(), "");
        bool areParamsEqual = true;
        for(int i = 0; i < prevCosmoParams_.size(); ++i)
        {
            if(prevCosmoParams_[i] != currentCosmoParams_[i])
            {
                areParamsEqual = false;
                break;
            }
        }

        needToCalculate = !areParamsEqual;
    }
    else
        prevCosmoParamsName_ = params.name();

    if(!needToCalculate)
        return;

    haveLow_ = false;
    haveLens_ = false;

    params_ = &params;
    prevCosmoParams_.swap(currentCosmoParams_);

    const bool wantT = true;
    const bool wantPol = (lowP_ || highP_ || lensingP_);
    const bool wantLens = lensingT_;

    cmb_->initialize(params, wantT, wantPol, true);

    std::vector<double>* tt = &clTT_;
    std::vector<double>* ee = (wantPol ? &clEE_ : NULL);
    std::vector<double>* te = (wantPol ? &clTE_ : NULL);
    std::vector<double>* bb = (lowP_ ? &clBB_ : NULL);

    cmb_->getLensedCl(tt, ee, te, bb);

    if(wantLens)
        cmb_->getCl(NULL, NULL, NULL, &clPP_, NULL, NULL);
}

void
PlanckLikelihood::setCls(const std::vector<double>* tt, const std::vector<double>* ee, const std::vector<double>* te, const std::vector<double> *bb, const std::vector<double>* pp)
{
    check(tt, "TT must be specified");

    clTT_ = *tt;

    if(ee)
        clEE_ = *ee;
    else
        clEE_.clear();

    if(te)
        clTE_ = *te;
    else
        clTE_.clear();

    if(bb)
        clBB_ = *bb;
    else
        clBB_.clear();

    if(pp)
        clPP_ = *pp;
    else
        clPP_.clear();

    haveLow_ = false;
    haveLens_ = false;
}

void
PlanckLikelihood::setAPlanck(double aPlanck)
{
    check(aPlanck > 0, "");
    aPlanck_ = aPlanck;

    haveLow_ = false;
    haveLens_ = false;
}

void
PlanckLikelihood::setAPol(double aPol)
{
    check(lowP_ || highP_ || lensingP_, "polarization not initialized");
    check(aPol > 0, "");
    aPol_ = aPol;

    haveLow_ = false;
    haveLens_ = false;
}

void
PlanckLikelihood::setHighExtraParams(const std::vector<double>& params)
{
    check(highT_, "high l likelihood not initialized");
    check(params.size() == (highP_ ? 32 : 15), "");
    check(highExtra_.size() == params.size(), "");
    for(int i = 0; i < highExtra_.size(); ++i)
        highExtra_[i] = params[i];
}

void
PlanckLikelihood::setBeamLeakageParams(const std::vector<double>& params)
{
    check(highP_, "high l polarization likelihood not initialized");
    check(params.size() == 60, "");
    check(beamExtra_.size() == 60, "");
    for(int i = 0; i < 60; ++i)
        beamExtra_[i] = params[i];
}

void
PlanckLikelihood::setSZPrior(bool szPrior)
{
    check(high_, "high likelihood not initialized");
    check(!highLikeLite_, "using the lite version of high likelihood");

    szPrior_ = szPrior;
}

double
PlanckLikelihood::lowLike()
{
    check(low_, "low l likelihood not initialized");

    if(haveLow_)
        return prevLow_;

    check(!clTT_.empty(), "Cl-s not computed");
    check(clTT_.size() >= lowLMax_ + 1, "");
    check(!lowP_ || !clTE_.empty(), "Cl-s not computed");
    check(!lowP_ || clTE_.size() >= lowLMax_ + 1, "");
    check(!lowP_ || !clEE_.empty(), "Cl-s not computed");
    check(!lowP_ || clEE_.size() >= lowLMax_ + 1, "");
    check(!lowP_ || !clBB_.empty(), "Cl-s not computed");
    check(!lowP_ || clBB_.size() >= lowLMax_ + 1, "");
    check(input_.size() >= 4 * (lowLMax_ + 1) + 1, "");
    output_screen2("Calculating low-l likelihood..." << std::endl);

    //Timer timer("PLANCK LOW-L LIKELIHOOD");
    //timer.start();

    std::vector<double>::iterator it = input_.begin();
    for(int l = 0; l <= lowLMax_; ++l)
    {
        *it = clTT_[l];
        ++it;
    }
    if(lowP_)
    {
        for(int l = 0; l <= lowLMax_; ++l)
        {
            *it = clEE_[l];
            ++it;
        }
        for(int l = 0; l <= lowLMax_; ++l)
        {
            *it = clBB_[l];
            ++it;
        }
        for(int l = 0; l <= lowLMax_; ++l)
        {
            *it = clTE_[l];
            ++it;
        }
    }
    *it = aPlanck_;
    const double l = clik_compute(low_, &(input_[0]), NULL);

    //timer.end();

    output_screen2("OK" << std::endl);

    prevLow_ = -2.0 * l;
    haveLow_ = true;

    return prevLow_;
}

double
PlanckLikelihood::highLike()
{
    check(high_, "high l likelihood not initialized");

    check(!clTT_.empty(), "Cl-s not computed");
    check(clTT_.size() >= highLMax_ + 1, "");
    check(!highP_ || !clTE_.empty(), "Cl-s not computed");
    check(!highP_ || clTE_.size() >= highLMax_ + 1, "");
    check(!highP_ || !clEE_.empty(), "Cl-s not computed");
    check(!highP_ || clEE_.size() >= highLMax_ + 1, "");
    check(input_.size() >= 3 * (highLMax_ + 1) + 100, "");
    output_screen2("Calculating high-l likelihood..." << std::endl);

    //Timer timer("PLANCK HIGH-L LIKELIHOOD");
    //timer.start();

    std::vector<double>::iterator it = input_.begin();
    for(int l = 0; l <= highLMax_; ++l)
    {
        *it = clTT_[l];
        ++it;
    }
    if(highP_)
    {
        for(int l = 0; l <= highLMax_; ++l)
        {
            *it = clEE_[l];
            ++it;
        }
        for(int l = 0; l <= highLMax_; ++l)
        {
            *it = clTE_[l];
            ++it;
        }
    }
    if(!highLikeLite_)
    {
        if(highP_)
        {
            check(highExtra_.size() == 32, "");
            for(int i = 0; i < 27; ++i)
            {
                *it = highExtra_[i];
                ++it;
            }
            check(beamExtra_.size() == 60, "");
            for(int i  = 0; i < 60; ++i)
            {
                *it = beamExtra_[i];
                ++it;
            }
            for(int i = 27; i < 32; ++i)
            {
                *it = highExtra_[i];
                ++it;
            }
            *it = aPol_;
            ++it;
        }
        else
        {
            check(highExtra_.size() == 15, "");
            for(int i = 0; i < 15; ++i)
            {
                *it = highExtra_[i];
                ++it;
            }
        }
    }
    *it = aPlanck_;
    double l = clik_compute(high_, &(input_[0]), NULL);

    if(szPrior_)
    {
        check(!highLikeLite_, "");
        // using 1.6*A_sz + k_sz = 9.5 +/- 3.0
        const double sz = 1.6 * highExtra_[3] + highExtra_[8];
        const double szMean = 9.5;
        const double szSigma = 3;
        l -= (sz - szMean) * (sz - szMean) / (2 * szSigma * szSigma);
    }

    //timer.end();
    output_screen2("OK" << std::endl);
    return -2.0 * l;
}

double
PlanckLikelihood::lensingLike()
{
    check(lens_, "lensing not initialized");

    if(haveLens_)
        return prevLens_;

    check(!clPP_.empty(), "Cl-s not computed");
    check(clPP_.size() >= lensLMax_ + 1, "");
    check(!clTT_.empty(), "Cl-s not computed");
    check(clTT_.size() >= lensLMax_ + 1, "");
    check(!lensingP_ || !clTE_.empty(), "Cl-s not computed");
    check(!lensingP_ || clTE_.size() >= lensLMax_ + 1, "");
    check(!lensingP_ || !clEE_.empty(), "Cl-s not computed");
    check(!lensingP_ || clEE_.size() >= lensLMax_ + 1, "");
    check(input_.size() >= 4 * (lensLMax_ + 1) + 1, "");
    output_screen2("Calculating lensing likelihood..." << std::endl);

    //Timer timer("PLANCK LENSING LIKELIHOOD");
    //timer.start();

    std::vector<double>::iterator it = input_.begin();
    for(int l = 0; l <= lensLMax_; ++l)
    {
        *it = clPP_[l];
        ++it;
    }
    for(int l = 0; l <= lensLMax_; ++l)
    {
        *it = clTT_[l];
        ++it;
    }
    if(lensingP_)
    {
        for(int l = 0; l <= lensLMax_; ++l)
        {
            *it = clEE_[l];
            ++it;
        }
        for(int l = 0; l <= lensLMax_; ++l)
        {
            *it = clTE_[l];
            ++it;
        }
    }
    *it = aPlanck_;
    const double l = clik_lensing_compute(static_cast<clik_lensing_object*>(lens_), &(input_[0]), NULL);

    //timer.end();
    output_screen2("OK" << std::endl);
    
    prevLens_ = -2.0 * l;
    haveLens_ = true;

    return prevLens_;
}


double
PlanckLikelihood::calculate(double* params, int nPar)
{
    //Timer timer("Planck likelihood timer");
    //timer.start();
    
    check(modelParams_, "model params must be set before calling this function");
    check(!vModel_.empty(), "");
    const int nModel = vModel_.size();
    
    check(nPar == (nModel + 1 + (highT_ && !highLikeLite_ ? (highP_ ? 33 : 15) : 0)), "");

    for(int i = 0; i < nModel; ++i)
        vModel_[i] = params[i];

    modelParams_->setAllParameters(vModel_);
    setCosmoParams(*modelParams_);

    aPlanck_ = params[nModel];

    if(highT_)
    {
        check(highExtra_.size() == (highP_ ? 32 : 15), "");
        for(int i = 0; i < highExtra_.size(); ++i)
            highExtra_[i] = params[nModel + 1 + i];
    }

    if(highP_)
        aPol_ = params[nModel + 33];

    //timer.end();
    return likelihood();
}

double
PlanckLikelihood::likelihood()
{
    double l = 0;
    if(low_) l += lowLike();
    if(high_) l += highLike();
    if(lens_) l += lensingLike();
    return l;
}

#else

namespace
{

struct PlanckLikelihoodContainer
{
    PlanckLikelihoodContainer() : commander(NULL), camspec(NULL), lens(NULL), pol(NULL), actspt(NULL)
    {
    }

    ~PlanckLikelihoodContainer()
    {
        if(commander)
        {
            clik_cleanup(&commander);
            commander = NULL;
        }
        if(camspec)
        {
            clik_cleanup(&camspec);
            camspec = NULL;
        }
        if(pol)
        {
            clik_cleanup(&pol);
            pol = NULL;
        }
        if(actspt)
        {
            clik_cleanup(&actspt);
            actspt = NULL;
        }
        if(lens)
        {
            clik_lensing_cleanup(&lens);
            lens = NULL;
        }
    }

    void* commander;
    void* camspec;
    void* lens;
    void* pol;
    void* actspt;
};

}

PlanckLikelihood::PlanckLikelihood(bool useCommander, bool useCamspec, bool useLensing, bool usePol, bool useActSpt, bool includeTensors, double kPerDecade, bool useOwnCmb) : spectraNames_(6), haveCommander_(false), havePol_(false), haveLens_(false), commander_(NULL), camspec_(NULL), lens_(NULL), pol_(NULL), actspt_(NULL), cmb_(NULL), modelParams_(NULL)
{
    check(useCommander || useCamspec || useLensing || usePol || useActSpt, "at least one likelihood must be specified");

    std::string planckLikeDir = PLANCK_DATA_DIR_STR;

    output_screen("Planck data dir = " << planckLikeDir << std::endl);

    int hasCl[6], lMax[6];
    parname* names;

    spectraNames_[0] = std::string("TT");
    spectraNames_[1] = std::string("EE");
    spectraNames_[2] = std::string("BB");
    spectraNames_[3] = std::string("TE");
    spectraNames_[4] = std::string("TB");
    spectraNames_[5] = std::string("EB");

    lMax_ = 0;
    commanderLMax_ = 0;
    camspecLMax_ = 0;
    polLMax_ = 0;
    lensingLMax_ = 0;
    actSptLMax_ = 0;

    static PlanckLikelihoodContainer cont;

    if(useCommander)
    {
        std::stringstream commanderPath;
        commanderPath << planckLikeDir << "/commander_v4.1_lm49.clik";
        char commanderPathCStr[100];
        std::strcpy(commanderPathCStr, commanderPath.str().c_str());

        if(!cont.commander)
            cont.commander = clik_init(commanderPathCStr, NULL);

        commander_ = cont.commander;

        clik_get_has_cl(commander_, hasCl, NULL);
        clik_get_lmax(commander_, lMax, NULL);
        for(int i = 0; i < 6; ++i)
        {
            output_screen1("Commander has " << spectraNames_[i] << " = " << hasCl[i] << " with l_max = " << lMax[i] << std::endl);
            if(lMax[i] > commanderLMax_)
                commanderLMax_ = lMax[i];
        }

        if(commanderLMax_ > lMax_)
            lMax_ = commanderLMax_;

        /*
        int extraParams = clik_get_extra_parameter_names(commander_, &names, NULL);
        output_screen1("Commander has " << extraParams << " extra parameters. Their names are as follows:" << std::endl);
        for(int i = 0; i < extraParams; ++i)
        {
            output_screen1(names[i] << std::endl);
        }
        output_screen1(std::endl);
        delete names;
        */
    }

    if(useCamspec)
    {
        std::stringstream camspecPath;
        camspecPath << planckLikeDir << "/CAMspec_v6.2TN_2013_02_26_dist.clik";
        char camspecPathCStr[100];
        std::strcpy(camspecPathCStr, camspecPath.str().c_str());

        if(!cont.camspec)
            cont.camspec = clik_init(camspecPathCStr, NULL);

        camspec_ = cont.camspec;

        clik_get_has_cl(camspec_, hasCl, NULL);
        clik_get_lmax(camspec_, lMax, NULL);
        for(int i = 0; i < 6; ++i)
        {
            output_screen1("Camspec has " << spectraNames_[i] << " = " << hasCl[i] << " with l_max = " << lMax[i] << std::endl);
            if(lMax[i] > camspecLMax_)
                camspecLMax_ = lMax[i];
        }

        if(camspecLMax_ > lMax_)
            lMax_ = camspecLMax_;

        /*
        int extraParams = clik_get_extra_parameter_names(camspec_, &names, NULL);
        output_screen1("Camspec has " << extraParams << " extra parameters. Their names are as follows:" << std::endl);
        for(int i = 0; i < extraParams; ++i)
        {
            output_screen1(names[i] << std::endl);
        }
        output_screen1(std::endl);
        delete names;
        */
    }

    if(usePol)
    {
        std::stringstream polPath;
        polPath << planckLikeDir << "/lowlike_v222.clik";
        char polPathCStr[100];
        std::strcpy(polPathCStr, polPath.str().c_str());

        if(!cont.pol)
            cont.pol = clik_init(polPathCStr, NULL);

        pol_ = cont.pol;

        clik_get_has_cl(pol_, hasCl, NULL);
        clik_get_lmax(pol_, lMax, NULL);
        for(int i = 0; i < 6; ++i)
        {
            output_screen1("Pol has " << spectraNames_[i] << " = " << hasCl[i] << " with l_max = " << lMax[i] << std::endl);
            if(lMax[i] > polLMax_)
                polLMax_ = lMax[i];
        }

        if(polLMax_ > lMax_)
            lMax_ = polLMax_;

        /*
        int extraParams = clik_get_extra_parameter_names(pol_, &names, NULL);
        output_screen1("Pol has " << extraParams << " extra parameters. Their names are as follows:" << std::endl);
        for(int i = 0; i < extraParams; ++i)
        {
            output_screen1(names[i] << std::endl);
        }
        output_screen1(std::endl);
        delete names;
        */
    }

    if(useLensing)
    {
        std::stringstream lensPath;
        lensPath << planckLikeDir << "/lensing_likelihood_v4_ref.clik_lensing";
        char lensPathCStr[100];
        std::strcpy(lensPathCStr, lensPath.str().c_str());

        if(!cont.lens)
            cont.lens = clik_lensing_init(lensPathCStr, NULL);

        lens_ = cont.lens;

        lensingLMax_ = clik_lensing_get_lmax(lens_,NULL);
        if(lensingLMax_ > lMax_)
            lMax_ = lensingLMax_;

        output_screen1("Lensing has l_max = " << lensingLMax_ << std::endl);

        /*
        int extraParams = clik_lensing_get_extra_parameter_names(lens_, &names, NULL);
        output_screen1("Lensing has " << extraParams << " extra parameters. Their names are as follows:" << std::endl);
        for(int i = 0; i < extraParams; ++i)
        {
            output_screen1(names[i] << std::endl);
        }
        output_screen1(std::endl);
        delete names;
        */
    }

    if(useActSpt)
    {
        std::stringstream actSptPath;
        actSptPath << planckLikeDir << "/actspt_2013_01.clik";
        char actSptPathCStr[100];
        std::strcpy(actSptPathCStr, actSptPath.str().c_str());

        if(!cont.actspt)
            cont.actspt = clik_init(actSptPathCStr, NULL);

        actspt_ = cont.actspt;

        clik_get_has_cl(actspt_, hasCl, NULL);
        clik_get_lmax(actspt_, lMax, NULL);
        for(int i = 0; i < 6; ++i)
        {
            output_screen1("ActSpt has " << spectraNames_[i] << " = " << hasCl[i] << " with l_max = " << lMax[i] << std::endl);
            if(lMax[i] > actSptLMax_)
                actSptLMax_ = lMax[i];
        }

        if(actSptLMax_ > lMax_)
            lMax_ = actSptLMax_;

        /*
        int extraParams = clik_get_extra_parameter_names(actspt_, &names, NULL);
        output_screen1("ActSpt has " << extraParams << " extra parameters. Their names are as follows:" << std::endl);
        for(int i = 0; i < extraParams; ++i)
        {
            output_screen1(names[i] << std::endl);
        }
        output_screen1(std::endl);
        delete names;
        */
    }

    output_screen1("Total l_max = " << lMax_ << std::endl);

    if(useOwnCmb)
    {
        cmb_ = new CMB;
        cmb_->preInitialize(lMax_ + 1000, false, true, includeTensors, lMax_ + 1000, kPerDecade);
    }
}

PlanckLikelihood::~PlanckLikelihood()
{
    if(cmb_)
        delete cmb_;
}

void
PlanckLikelihood::setCosmoParams(const CosmologicalParams& params)
{
    check(cmb_, "own CMB must be used (set in constructor)");

    bool needToCalculate = true;
    params.getAllParameters(currentCosmoParams_);
    if(params.name() == prevCosmoParamsName_)
    {
        check(currentCosmoParams_.size() == prevCosmoParams_.size(), "");
        bool areParamsEqual = true;
        for(int i = 0; i < prevCosmoParams_.size(); ++i)
        {
            if(prevCosmoParams_[i] != currentCosmoParams_[i])
            {
                areParamsEqual = false;
                break;
            }
        }

        needToCalculate = !areParamsEqual;
    }
    else
        prevCosmoParamsName_ = params.name();

    if(!needToCalculate)
        return;

    haveCommander_ = false;
    havePol_ = false;
    haveLens_ = false;

    params_ = &params;
    prevCosmoParams_.swap(currentCosmoParams_);

    bool wantT = true;
    bool wantPol = (pol_ != NULL);
    bool wantLens = (lens_ != NULL);

    cmb_->initialize(params, wantT, wantPol, true);

    std::vector<double>* tt = &clTT_;
    std::vector<double>* ee = (wantPol ? &clEE_ : NULL);
    std::vector<double>* te = (wantPol ? &clTE_ : NULL);

    cmb_->getLensedCl(tt, ee, te);

    if(wantLens)
        cmb_->getCl(NULL, NULL, NULL, &clPP_, NULL, NULL);
}

void
PlanckLikelihood::setCls(const std::vector<double>* tt, const std::vector<double>* ee, const std::vector<double>* te, const std::vector<double>* pp)
{
    check(tt, "TT must be specified");
    //check(ee || !pol_, "EE must be specified since polarization likelihood is used");
    //check(te || !pol_, "TE must be specified since polarization likelihood is used");
    //check(pp || !lens_, "PP must be specified since lensing likelihood is used");

    clTT_ = *tt;

    if(ee)
        clEE_ = *ee;
    else
        clEE_.clear();

    if(te)
        clTE_ = *te;
    else
        clTE_.clear();

    if(pp)
        clPP_ = *pp;
    else
        clPP_.clear();

    haveCommander_ = false;
    havePol_ = false;
    haveLens_ = false;
}

void
PlanckLikelihood::setCamspecExtraParams(double A_ps_100, double A_ps_143, double A_ps_217, double A_cib_143, double A_cib_217, double A_sz, double r_ps, double r_cib, double n_Dl_cib, double cal_100, double cal_217, double xi_sz_cib, double A_ksz, double Bm_1_1)
{
    check(camspec_, "Camspec likelihood not initialized");

    camspecExtra_.clear();
    camspecExtra_.push_back(A_ps_100);
    camspecExtra_.push_back(A_ps_143);
    camspecExtra_.push_back(A_ps_217);
    camspecExtra_.push_back(A_cib_143);
    camspecExtra_.push_back(A_cib_217);
    camspecExtra_.push_back(A_sz);
    camspecExtra_.push_back(r_ps);
    camspecExtra_.push_back(r_cib);
    camspecExtra_.push_back(n_Dl_cib);
    camspecExtra_.push_back(cal_100);
    camspecExtra_.push_back(cal_217);
    camspecExtra_.push_back(xi_sz_cib);
    camspecExtra_.push_back(A_ksz);
    camspecExtra_.push_back(Bm_1_1);
}

void
PlanckLikelihood::setActSptExtraParams(double A_sz, double A_ksz, double xi_sz_cib, double a_ps_act_148, double a_ps_act_217, double a_ps_spt_95, double a_ps_spt_150, double a_ps_spt_220, double A_cib_143, double A_cib_217, double n_Dl_cib, double r_ps_spt_95x150, double r_ps_spt_95x220, double r_ps_150x220, double r_cib, double a_gs, double a_ge, double cal_acts_148, double cal_acts_217, double cal_acte_148, double cal_acte_217, double cal_spt_95, double cal_spt_150, double cal_spt_220)
{
    check(actspt_, "ActSpt likelihood not initialized");

    actSptExtra_.clear();
    actSptExtra_.push_back(A_sz);
    actSptExtra_.push_back(A_ksz);
    actSptExtra_.push_back(xi_sz_cib);
    actSptExtra_.push_back(a_ps_act_148);
    actSptExtra_.push_back(a_ps_act_217);
    actSptExtra_.push_back(a_ps_spt_95);
    actSptExtra_.push_back(a_ps_spt_150);
    actSptExtra_.push_back(a_ps_spt_220);
    actSptExtra_.push_back(A_cib_143);
    actSptExtra_.push_back(A_cib_217);
    actSptExtra_.push_back(n_Dl_cib);
    actSptExtra_.push_back(r_ps_spt_95x150);
    actSptExtra_.push_back(r_ps_spt_95x220);
    actSptExtra_.push_back(r_ps_150x220);
    actSptExtra_.push_back(r_cib);
    actSptExtra_.push_back(a_gs);
    actSptExtra_.push_back(a_ge);
    actSptExtra_.push_back(cal_acts_148);
    actSptExtra_.push_back(cal_acts_217);
    actSptExtra_.push_back(cal_acte_148);
    actSptExtra_.push_back(cal_acte_217);
    actSptExtra_.push_back(cal_spt_95);
    actSptExtra_.push_back(cal_spt_150);
    actSptExtra_.push_back(cal_spt_220);
}

double
PlanckLikelihood::commanderLike()
{
    check(commander_, "commander not initialized");

    if(haveCommander_)
        return prevCommander_;

    check(!clTT_.empty(), "Cl-s not computed");

    check(clTT_.size() >= commanderLMax_ + 1, "");

    output_screen2("Calculating commander likelihood..." << std::endl);
    const double l = clik_compute(commander_, &(clTT_[0]), NULL);
    output_screen2("OK" << std::endl);

    prevCommander_ = -2.0 * l;
    haveCommander_ = true;

    return -2.0 * l;
}

double
PlanckLikelihood::camspecLike()
{
    check(camspec_, "camspec not initialized");
    check(!clTT_.empty(), "Cl-s not computed");
    check(!camspecExtra_.empty(), "camspec extra parameters not initialized");

    check(clTT_.size() >= camspecLMax_ + 1, "");
    check(camspecExtra_.size() == 14, "");

    output_screen2("Calculating camspec likelihood..." << std::endl);
    std::vector<double> input;
    input.insert(input.end(), clTT_.begin(), clTT_.begin() + camspecLMax_ + 1);
    input.insert(input.end(), camspecExtra_.begin(), camspecExtra_.end());

    const double l = clik_compute(camspec_, &(input[0]), NULL);
    output_screen2("OK" << std::endl);
    return -2.0 * l;
}

double
PlanckLikelihood::polLike()
{
    check(pol_, "pol not initialized");
    
    if(havePol_)
        return prevPol_;

    check(!clTT_.empty(), "Cl-s not computed");
    check(!clEE_.empty(), "EE Cl-s missing");
    check(!clTE_.empty(), "TE Cl-s missing");

    check(clTT_.size() >= polLMax_ + 1, "");
    check(clEE_.size() >= polLMax_ + 1, "");
    check(clTE_.size() >= polLMax_ + 1, "");

    output_screen2("Calculating polarization likelihood..." << std::endl);
    std::vector<double> input;
    input.insert(input.end(), clTT_.begin(), clTT_.begin() + polLMax_ + 1);
    input.insert(input.end(), clEE_.begin(), clEE_.begin() + polLMax_ + 1);
    input.insert(input.end(), polLMax_ + 1, 0.0);
    input.insert(input.end(), clTE_.begin(), clTE_.begin() + polLMax_ + 1);

    const double l = clik_compute(pol_, &(input[0]), NULL);
    output_screen2("OK" << std::endl);

    prevPol_ = -2.0 * l;
    havePol_ = true;

    return -2.0 * l;
}

double
PlanckLikelihood::lensingLike()
{
    check(lens_, "lensing not initialized");

    if(haveLens_)
        return prevLens_;

    check(!clTT_.empty(), "Cl-s not computed");
    check(!clPP_.empty(), "Lensing Cl-s missing");

    check(clTT_.size() >= lensingLMax_ + 1, "");
    check(clPP_.size() >= lensingLMax_ + 1, "");

    output_screen2("Calculating lensing likelihood..." << std::endl);
    std::vector<double> input;
    input.insert(input.end(), clPP_.begin(), clPP_.begin() + lensingLMax_ + 1);
    input.insert(input.end(), clTT_.begin(), clTT_.begin() + lensingLMax_ + 1);
    
    const double l = clik_lensing_compute(lens_, &(input[0]), NULL);
    output_screen2("OK" << std::endl);

    prevLens_ = -2.0 * l;
    haveLens_ = true;

    return -2.0 * l;
}

double
PlanckLikelihood::actSptLike()
{
    check(actspt_, "actspt not initialized");
    check(!clTT_.empty(), "Cl-s not computed");
    check(!actSptExtra_.empty(), "actspt extra parameters not initialized");

    check(clTT_.size() >= actSptLMax_ + 1, "");
    check(actSptExtra_.size() == 24, "");

    output_screen2("Calculating act-spt likelihood..." << std::endl);
    std::vector<double> input;
    input.insert(input.end(), clTT_.begin(), clTT_.begin() + actSptLMax_ + 1);
    input.insert(input.end(), actSptExtra_.begin(), actSptExtra_.end());

    const double l = clik_compute(actspt_, &(input[0]), NULL);
    output_screen2("OK" << std::endl);

    return -2.0 * l;
}

double
PlanckLikelihood::calculate(double* params, int nPar)
{
    //Timer timer("Planck likelihood timer");
    //timer.start();
    
    check(modelParams_, "model params must be set before calling this function");
    check(!vModel_.empty(), "");
    const int nModel = vModel_.size();
    
    check(nPar == (nModel + (camspec_ ? 14 : 0) + (actspt_ ? 24 : 0)), "");

    for(int i = 0; i < nModel; ++i)
        vModel_[i] = params[i];

    modelParams_->setAllParameters(vModel_);
    setCosmoParams(*modelParams_);

    if(camspec_)
    {
        camspecExtra_.clear();
        camspecExtra_.insert(camspecExtra_.end(), params + nModel, params + nModel + 14);
    }

    if(actspt_)
    {
        int paramShift = nModel + (camspec_ ? 14 : 0);
        actSptExtra_.clear();
        actSptExtra_.insert(actSptExtra_.end(), params + paramShift, params + paramShift + 24);
    }

    //timer.end();
    return likelihood();
}

double
PlanckLikelihood::likelihood()
{
    double l = 0;

    if(commander_)
        l += commanderLike();

    if(camspec_)
        l += camspecLike();

    if(pol_)
        l += polLike();

    if(lens_)
        l += lensingLike();

    if(actspt_)
        l+= actSptLike();

    return l;
}

#endif

