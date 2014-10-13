#ifdef COSMO_MPI
#include <mpi.h>
#endif

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

PlanckLikelihood::PlanckLikelihood(bool useCommander, bool useCamspec, bool useLensing, bool usePol, bool useActSpt, bool includeTensors, double kPerDecade) : spectraNames_(6), cosmoParams_(6), prevCosmoCalculated_(false), needToCalculate_(true), commander_(NULL), camspec_(NULL), lens_(NULL), pol_(NULL), actspt_(NULL)
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
        if(!cont.commander) cont.commander = clik_init(commanderPathCStr, NULL);
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
        if(!cont.camspec) cont.camspec = clik_init(camspecPathCStr, NULL);
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
        if(!cont.pol) cont.pol = clik_init(polPathCStr, NULL);
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
        if(!cont.lens) cont.lens = clik_lensing_init(lensPathCStr, NULL);
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
        if(!cont.actspt) cont.actspt = clik_init(actSptPathCStr, NULL);
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
    cmb_.preInitialize(lMax_ + 1000, false, true, includeTensors, lMax_ + 1000, kPerDecade);

    useCMB_ = &cmb_;
}

void
PlanckLikelihood::setCosmoParams(const CosmologicalParams& params)
{
    params_ = &params;

    bool wantT = true;
    bool wantPol = (pol_ != NULL);
    bool wantLens = (lens_ != NULL);

    cmb_.initialize(params, wantT, wantPol, true);
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

void
PlanckLikelihood::calculateCls()
{
    bool wantT = true;
    bool wantPol = (pol_ != NULL);
    bool wantLens = (lens_ != NULL);

    std::vector<double>* tt = &clTT_;
    std::vector<double>* ee = (wantPol ? &clEE_ : NULL);
    std::vector<double>* te = (wantPol ? &clTE_ : NULL);
    
    useCMB_->getLensedCl(tt, ee, te);

    if(wantLens)
        useCMB_->getCl(NULL, NULL, NULL, &clPP_, NULL, NULL);

    // to be removed
    /*
    std::ifstream in("l_list.txt");
    if(!in)
    {
        StandardException exc;
        std::string exceptionStr = "Cannot open input file l_list.txt.";
        exc.set(exceptionStr);
        throw exc;
    }

    int totalNum;
    in >> totalNum;
    std::vector<int> lList(totalNum);
    for(int i = 0; i < totalNum; ++i)
        in >> lList[i];
    in.close();

    std::stringstream nameStr;
    nameStr << "cl_data";
#ifdef COSMO_MPI
    int mpif;
    MPI_Initialized(&mpif);
    if(mpif)
    {
        int n;
        MPI_Comm_size(MPI_COMM_WORLD, &n);
        if(n > 1)
        {
            int p;
            MPI_Comm_rank(MPI_COMM_WORLD, &p);

            nameStr << '_' << p;
        }
    }
#endif
    nameStr << ".dat";

    std::vector<double> buff;
    buff.push_back(params_->getOmBH2());
    buff.push_back(params_->getOmCH2());
    buff.push_back(params_->getH());
    buff.push_back(params_->getTau());
    buff.push_back(params_->getNs());
    buff.push_back(params_->getAs());

    for(int i = 0; i < totalNum; ++i)
    {
        const int l = lList[i];
        buff.push_back(clTT_[l]);
        buff.push_back(clEE_[l]);
        buff.push_back(clTE_[l]);
        buff.push_back(clPP_[l]);
    }

    std::ofstream out(nameStr.str().c_str(), std::ios::binary | std::ios::app);
    out.write((char*)(&(buff[0])), buff.size() * sizeof(double));
    out.close();
    */
}

double
PlanckLikelihood::commanderLike()
{
    check(commander_, "commander not initialized");

    if(!needToCalculate_)
        return prevCommander_;

    check(!clTT_.empty(), "Cl-s not computed");

    check(clTT_.size() >= commanderLMax_ + 1, "");

    output_screen2("Calculating commander likelihood..." << std::endl);
    const double l = clik_compute(commander_, &(clTT_[0]), NULL);
    output_screen2("OK" << std::endl);
    prevCommander_ = -2.0 * l;
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
    
    if(!needToCalculate_)
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

    return -2.0 * l;
}

double
PlanckLikelihood::lensingLike()
{
    check(lens_, "lensing not initialized");

    if(!needToCalculate_)
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
    
    check(useCMB_ == &cmb_, "");

    check(nPar == (6 + (camspec_ ? 14 : 0) + (actspt_ ? 24 : 0)), "");
    const double pivot = 0.05;

    needToCalculate_ = true;
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

        needToCalculate_ = !areParamsEqual;
    }

    if(needToCalculate_)
    {
        LambdaCDMParams lcdmParams(params[0], params[1], params[2], params[3], params[4], std::exp(params[5]) / 1e10, pivot);
        setCosmoParams(lcdmParams);
    }

    if(camspec_)
    {
        camspecExtra_.clear();
        camspecExtra_.insert(camspecExtra_.end(), params + 6, params + 20);
    }

    if(actspt_)
    {
        int paramShift = (camspec_ ? 20 : 6);
        actSptExtra_.clear();
        actSptExtra_.insert(actSptExtra_.end(), params + paramShift, params + paramShift + 24);
    }

    if(needToCalculate_)
    {
        calculateCls();
        check(cosmoParams_.size() == 6, "");
        for(int i = 0; i < 6; ++i)
            cosmoParams_[i] = params[i];
        prevCosmoCalculated_ = true;
    }

    //timer.end();
    return likelihood();
}

double
PlanckLikelihood::likelihood()
{
    double l = 0;
    if(commander_) l += commanderLike();
    if(camspec_) l += camspecLike();
    if(pol_) l += polLike();
    if(lens_) l += lensingLike();
    if(actspt_) l+= actSptLike();

    return l;
}

PlanckLikelihood::~PlanckLikelihood()
{
}

