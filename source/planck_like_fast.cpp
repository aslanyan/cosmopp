#include <cosmo_mpi.hpp>

#include <string>
#include <sstream>
#include <iomanip>
#include <ctime>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <cubic_spline.hpp>
#include <planck_like_fast.hpp>

namespace
{

#ifdef COSMO_PLANCK_15

class PlanckLikeFastClFunc : public Math::RealFunctionMultiToMulti
{
public:
    PlanckLikeFastClFunc(CosmologicalParams* params, CMB* cmb, const std::vector<double>& lList, bool usePolarization, bool useBB, bool useLensing) : cosmoParams_(params), cmb_(cmb), lList_(lList), usePol_(usePolarization), useBB_(useBB), useLensing_(useLensing), clTT_(NULL), clEE_(NULL), clTE_(NULL), clBB_(NULL), clPP_(NULL)
    {
        check(!lList_.empty(), "");

        clTT_ = new std::vector<double>;

        if(usePol_)
        {
            clEE_ = new std::vector<double>;
            clTE_ = new std::vector<double>;
        }

        if(useBB_)
            clBB_ = new std::vector<double>;

        if(useLensing_)
            clPP_ = new std::vector<double>;
    }

    ~PlanckLikeFastClFunc()
    {
        if(clTT_) delete clTT_;
        if(clEE_) delete clEE_;
        if(clTE_) delete clTE_;
        if(clBB_) delete clBB_;
        if(clPP_) delete clPP_;
    }

    virtual void evaluate(const std::vector<double>& x, std::vector<double>* res) const
    {
#ifdef CHECKS_ON
        const bool success =
#endif
        cosmoParams_->setAllParameters(x);
        
        check(success, "");

        cmb_->initialize(*cosmoParams_, true, usePol_, true);

        cmb_->getLensedCl(clTT_, clEE_, clTE_, clBB_);
        if(useLensing_)
            cmb_->getCl(NULL, NULL, NULL, clPP_);

        int nCl = 1;
        if(usePol_) nCl += 2;
        if(useBB_) ++nCl;
        if(useLensing_) ++nCl;

        res->resize(nCl * lList_.size());
        std::vector<double>::iterator it = res->begin();
        for(int i = 0; i < lList_.size(); ++i)
        {
            check(it < res->end(), "");
            *(it++) = (*clTT_)[lList_[i]];
        }

        if(usePol_)
        {
            for(int i = 0; i < lList_.size(); ++i)
            {
                check(it < res->end(), "");
                *(it++) = (*clEE_)[lList_[i]];
            }
            for(int i = 0; i < lList_.size(); ++i)
            {
                check(it < res->end(), "");
                *(it++) = (*clTE_)[lList_[i]];
            }
        }

        if(useBB_)
        {
            for(int i = 0; i < lList_.size(); ++i)
            {
                check(it < res->end(), "");
                *(it++) = (*clBB_)[lList_[i]];
            }
        }

        if(useLensing_)
        {
            for(int i = 0; i < lList_.size(); ++i)
            {
                check(it < res->end(), "");
                *(it++) = (*clPP_)[lList_[i]];
            }
        }

        check(it == res->end(), "");
    }

private:
    CosmologicalParams* cosmoParams_;
    CMB* cmb_;
    const std::vector<double>& lList_;
    bool usePol_, useBB_, useLensing_;

    std::vector<double> *clTT_, *clEE_, *clTE_, *clBB_, *clPP_;
};

class PlanckLikeFastClErrorFunc : public Math::RealFunctionMultiDim
{
public:
    PlanckLikeFastClErrorFunc(PlanckLikelihood* like, const std::vector<double>& lList, bool usePolarization, bool useBB, bool useLensing, bool highP) : like_(like), lList_(lList), usePol_(usePolarization), useBB_(useBB), useLensing_(useLensing)
    {
        check(!lList_.empty(), "");
        tt_.resize(lList_.size(), 0);
        clTT_.resize(like_->getLMax() + 1, 0);

        if(usePol_)
        {
            ee_.resize(lList_.size(), 0);
            te_.resize(lList_.size(), 0);

            clEE_.resize(like_->getLMax() + 1, 0);
            clTE_.resize(like_->getLMax() + 1, 0);
        }

        if(useBB_)
        {
            bb_.resize(lList_.size(), 0);
            clBB_.resize(like_->getLMax() + 1, 0);
        }

        if(useLensing_)
        {
            pp_.resize(lList_.size(), 0);
            clPP_.resize(like_->getLMax() + 1, 0);
        }

        extraParams_.resize((highP ? 32 : 15));
        extraParams_[0] = 100;
        extraParams_[1] = -1.3;
        extraParams_[2] = 0.5;
        extraParams_[3] = 5;
        extraParams_[4] = 200;
        extraParams_[5] = 200;
        extraParams_[6] = 200;
        extraParams_[7] = 200;
        extraParams_[8] = 5;
        extraParams_[9] = 7;
        extraParams_[10] = 9;
        extraParams_[11] = 21;
        extraParams_[12] = 80;

        if(highP)
        {
            extraParams_[13] = 0.06;
            extraParams_[14] = 0.05;
            extraParams_[15] = 0.11;
            extraParams_[16] = 0.10;
            extraParams_[17] = 0.24;
            extraParams_[18] = 0.72;
            extraParams_[19] = -2.4;
            extraParams_[20] = 0.14;
            extraParams_[21] = 0.12;
            extraParams_[22] = 0.30;
            extraParams_[23] = 0.24;
            extraParams_[24] = 0.60;
            extraParams_[25] = 1.8;
            extraParams_[26] = -2.4;
            extraParams_[27] = 0.999;
            extraParams_[28] = 0.995;
            extraParams_[29] = 1;
            extraParams_[30] = 1;
            extraParams_[31] = 1;
        }
        else
        {
            extraParams_[13] = 0.999;
            extraParams_[14] = 0.995;
        }
    }

    ~PlanckLikeFastClErrorFunc() {}

    virtual double evaluate(const std::vector<double>& x) const
    {
        int nCl = 1;
        if(usePol_) nCl += 2;
        if(useBB_) ++nCl;
        if(useLensing_) ++nCl;

        check(x.size() == nCl * lList_.size(), "");

        std::vector<double>::const_iterator it = x.begin();

        check(tt_.size() == lList_.size(), "");
        for(int i = 0; i < tt_.size(); ++i)
        {
            check(it < x.end(), "");
            const double l = lList_[i];
            (*const_cast<std::vector<double>* >(&tt_))[i] = *(it++) * l * (l + 1);
        }

        if(usePol_)
        {
            check(ee_.size() == lList_.size(), "");
            check(te_.size() == lList_.size(), "");

            for(int i = 0; i < ee_.size(); ++i)
            {
                check(it < x.end(), "");
                const double l = lList_[i];
                (*const_cast<std::vector<double>* >(&ee_))[i] = *(it++) * l * (l + 1);
            }

            for(int i = 0; i < te_.size(); ++i)
            {
                check(it < x.end(), "");
                const double l = lList_[i];
                (*const_cast<std::vector<double>* >(&te_))[i] = *(it++) * l * (l + 1);
            }
        }

        if(useBB_)
        {
            check(bb_.size() == lList_.size(), "");

            for(int i = 0; i < bb_.size(); ++i)
            {
                check(it < x.end(), "");
                const double l = lList_[i];
                (*const_cast<std::vector<double>* >(&bb_))[i] = *(it++) * l * (l + 1);
            }
        }

        if(useLensing_)
        {
            check(pp_.size() == lList_.size(), "");

            for(int i = 0; i < pp_.size(); ++i)
            {
                check(it < x.end(), "");
                const double l = lList_[i];
                (*const_cast<std::vector<double>* >(&pp_))[i] = *(it++) * l * l * l * l;
            }
        }

        check(it == x.end(), "");

        Math::CubicSpline cs(lList_, tt_);
        const int lMax = like_->getLMax();
        check(clTT_.size() == lMax + 1, "");
        for(int l = 2; l <= lMax; ++l)
            (*const_cast<std::vector<double>* >(&clTT_))[l] = cs.evaluate(double(l)) / (l * (l + 1));

        if(usePol_)
        {
            Math::CubicSpline cs1(lList_, ee_);
            check(clEE_.size() == lMax + 1, "");
            for(int l = 2; l <= lMax; ++l)
                (*const_cast<std::vector<double>* >(&clEE_))[l] = cs1.evaluate(double(l)) / (l * (l + 1));

            Math::CubicSpline cs2(lList_, ee_);
            check(clTE_.size() == lMax + 1, "");
            for(int l = 2; l <= lMax; ++l)
                (*const_cast<std::vector<double>* >(&clTE_))[l] = cs2.evaluate(double(l)) / (l * (l + 1));
        }

        if(useBB_)
        {
            Math::CubicSpline cs1(lList_, bb_);
            check(clBB_.size() == lMax + 1, "");
            for(int l = 2; l <= lMax; ++l)
                (*const_cast<std::vector<double>* >(&clBB_))[l] = cs1.evaluate(double(l)) / (l * (l + 1));
        }

        if(useLensing_)
        {
            Math::CubicSpline cs1(lList_, pp_);
            check(clPP_.size() == lMax + 1, "");
            for(int l = 2; l <= lMax; ++l)
                (*const_cast<std::vector<double>* >(&clPP_))[l] = cs1.evaluate(double(l)) / (l * l * l * l);
        }

        like_->setCls(&clTT_, (usePol_ ? &clEE_ : NULL), (usePol_ ? &clTE_ : NULL), (useBB_ ? &clBB_ : NULL), (useLensing_ ? &clPP_ : NULL)); 
        like_->setHighExtraParams(extraParams_);
        like_->setAPlanck(1);
        like_->setAPol(1);

        return like_->likelihood();
    }

private:
    PlanckLikelihood* like_;
    const std::vector<double>& lList_;
    const bool usePol_, useBB_, useLensing_;

    std::vector<double> clTT_, clEE_, clTE_, clBB_, clPP_;
    std::vector<double> tt_, ee_, te_, bb_, pp_;

    std::vector<double> extraParams_;
};

class PlanckLikeFastFunc : public Math::RealFunctionMultiToMulti
{
public:
    PlanckLikeFastFunc(PlanckLikelihood *like) : like_(like)
    {
    }

    ~PlanckLikeFastFunc()
    {
    }

    virtual void evaluate(const std::vector<double>& x, std::vector<double>* res) const
    {
        res->resize(1, like_->calculate(const_cast<double*>(&(x[0])), x.size()));
    }

private:
    PlanckLikelihood* like_;
};

class PlanckLikeFastErrorFunc : public Math::RealFunctionMultiDim
{
public:
    PlanckLikeFastErrorFunc() {}

    virtual double evaluate(const std::vector<double>& x) const
    {
        check(x.size() == 1, "");
        return x[0];
    }
};

#else

class PlanckLikeFastFunc : public Math::RealFunctionMultiToMulti
{
public:
    PlanckLikeFastFunc(CosmologicalParams* params, PlanckLikelihood* like, CMB* cmb, const std::vector<double>& lList, bool useCommander, bool useLensing, bool usePolarization) : cosmoParams_(params), like_(like), cmb_(cmb), lList_(lList), useCommander_(useCommander), useLensing_(useLensing), usePol_(usePolarization), clTT_(NULL), clEE_(NULL), clTE_(NULL), clPP_(NULL)
    {
        check(!lList_.empty(), "");

        clTT_ = new std::vector<double>;
        if(useLensing_)
            clPP_ = new std::vector<double>;

        if(usePol_)
        {
            clEE_ = new std::vector<double>;
            clTE_ = new std::vector<double>;
        }
    }

    ~PlanckLikeFastFunc()
    {
        if(clTT_) delete clTT_;
        if(clEE_) delete clEE_;
        if(clTE_) delete clTE_;
        if(clPP_) delete clPP_;
    }

    virtual void evaluate(const std::vector<double>& x, std::vector<double>* res) const
    {
#ifdef CHECKS_ON
        const bool success =
#endif
        cosmoParams_->setAllParameters(x);
        
        check(success, "");

        cmb_->initialize(*cosmoParams_, true, usePol_, true);

        cmb_->getLensedCl(clTT_, clEE_, clTE_);
        if(useLensing_)
            cmb_->getCl(NULL, NULL, NULL, clPP_);

        like_->setCls(clTT_, clEE_, clTE_, clPP_);

        res->resize(3 + lList_.size());
        for(int i = 0; i < lList_.size(); ++i)
            (*res)[i] = (*clTT_)[lList_[i]];

        (*res)[lList_.size()] = (useCommander_ ? like_->commanderLike() : 0);
        (*res)[lList_.size() + 1] = (useLensing_ ? like_->lensingLike() : 0);
        (*res)[lList_.size() + 2] = (usePol_ ? like_->polLike() : 0);
    }

private:
    CosmologicalParams* cosmoParams_;
    PlanckLikelihood* like_;
    CMB* cmb_;
    const std::vector<double>& lList_;
    bool useCommander_;
    bool useLensing_;
    bool usePol_;

    std::vector<double> *clTT_, *clEE_, *clTE_, *clPP_;
};

class PlanckLikeFastErrorFunc : public Math::RealFunctionMultiDim
{
public:
    PlanckLikeFastErrorFunc(PlanckLikelihood* like, const std::vector<double>& lList, bool useCamspec, bool useActSpt) : like_(like), lList_(lList), useCamspec_(useCamspec), useActSpt_(useActSpt)
    {
        check(!lList_.empty(), "");
        cl_.resize(lList_.size(), 0);
        clTT_.resize(like_->getLMax() + 1, 0);
    }

    ~PlanckLikeFastErrorFunc() {}

    virtual double evaluate(const std::vector<double>& x) const
    {
        check(x.size() == lList_.size() + 3, "");

        double res = x[lList_.size()] + x[lList_.size() + 1] + x[lList_.size() + 2];

        if(!useCamspec_ && !useActSpt_)
            return res;

        check(cl_.size() == lList_.size(), "");
        for(int i = 0; i < cl_.size(); ++i)
        {
            const double l = lList_[i];
            (*const_cast<std::vector<double>* >(&cl_))[i] = x[i] * l * (l + 1);
        }

        Math::CubicSpline cs(lList_, cl_);

        const int lMax = like_->getLMax();
        check(clTT_.size() == lMax + 1, "");
        for(int l = 2; l <= lMax; ++l)
            (*const_cast<std::vector<double>* >(&clTT_))[l] = cs.evaluate(double(l)) / (l * (l + 1));

        like_->setCls(&clTT_);

        if(useCamspec_)
        {
            //like_->setCamspecExtraParams(100, 50, 100, 10, 30, 5, 0.9, 0.4, 0.5, 1, 1, 0.5, 5, 0.5);
            like_->setCamspecExtraParams(168, 66, 98, 3.5, 34.5, 4.9, 0.971, 0.37, 0.65, 1.0005, 0.996, 0.47, 1.5, 0.6);
            res += like_->camspecLike();
        }

        if(useActSpt_)
        {
            like_->setActSptExtraParams(5, 5, 0.5, 15, 100, 15, 15, 100, 10, 30, 0.5, 0.5, 0.5, 0.5, 0.4, 0.4, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
            res += like_->actSptLike();
        }

        return res;
    }

private:
    PlanckLikelihood* like_;
    const std::vector<double>& lList_;
    bool useCamspec_;
    bool useActSpt_;

    std::vector<double> cl_, clTT_;
};

#endif

} // namespace

#ifdef COSMO_PLANCK_15

PlanckLikeFast::PlanckLikeFast(CosmologicalParams* params, bool lowT, bool lowP, bool highT, bool highP, bool highLikeLite, bool lensingT, bool lensingP, bool includeTensors, double kPerDecade, double precision, unsigned long minCount) : cosmoParams_(params), lowT_(lowT), lowP_(lowP), highT_(highT), highP_(highP), highLikeLite_(highLikeLite), lensingT_(lensingT), lensingP_(lensingP), like_(lowT, lowP, highT, highP, highLikeLite, lensingT, lensingP, includeTensors, kPerDecade, !highT || highLikeLite), logError_(false), rand_(std::time(0), 0, 1)
{
    check(!lowP || lowT, "cannot include lowP without lowT");
    check(!highP || highT, "cannot include highP without highT");
    check(!lensingP || lensingT, "cannot include lensingP wihtout lensingT");

    check(lowT || highT || lensingT, "at least one likelihood must be specified");

    cosmoParams_->getAllParameters(cosmoParamsVec_);
    check(!cosmoParamsVec_.empty(), "");

    nParams_ = cosmoParamsVec_.size() + 1;

    lMax_ = like_.getLMax();

    if(highT_ && !highLikeLite_)
    {
        nParams_ += (highP_ ? 33 : 15);

        const bool usePol = (lowP || highP || lensingP);
        const bool useBB = lowP;
        const bool useLensing = lensingT;

        cmb_.preInitialize(lMax_ + 1000, false, true, includeTensors, lMax_ + 1000, kPerDecade);
        cmb_.initialize(*cosmoParams_, true, usePol, true);
        cmb_.getLList(lList_);

        std::vector<double>::iterator it = std::lower_bound(lList_.begin(), lList_.end(), double(lMax_));
        check(it != lList_.end(), "the l list is too short");
        ++it;
        if(it != lList_.end())
            lList_.erase(it, lList_.end());

        PlanckLikeFastClFunc* f = new PlanckLikeFastClFunc(cosmoParams_, &cmb_, lList_, usePol, useBB, useLensing);
        func_ = f;

        PlanckLikeFastClErrorFunc* errorFunc = new PlanckLikeFastClErrorFunc(&like_, lList_, usePol, useBB, useLensing, highP);
        errorFunc_ = errorFunc;

        std::stringstream fileName;
        fileName << "planck_fast_" << cosmoParams_->name() << ".dat";

        int nCl = 1;
        if(usePol) nCl += 2;
        if(useBB) ++nCl;
        if(useLensing) ++nCl;

        layg_ = new LearnAsYouGo(cosmoParamsVec_.size(), nCl * lList_.size(), *f, *errorFunc, minCount, precision, fileName.str().c_str());

        layg_->logIntoFile("planck_like_fast_log");

        tt_.resize(lList_.size(), 0);
        clTT_.resize(lMax_ + 1, 0);

        if(usePol)
        {
            ee_.resize(lList_.size(), 0);
            te_.resize(lList_.size(), 0);

            clEE_.resize(lMax_ + 1, 0);
            clTE_.resize(lMax_ + 1, 0);
        }

        if(useBB)
        {
            bb_.resize(lList_.size(), 0);
            clBB_.resize(lMax_ + 1, 0);
        }

        if(useLensing)
        {
            pp_.resize(lList_.size(), 0);
            clPP_.resize(lList_.size(), 0);
        }

        extraParams_.resize(highP_ ? 32 : 15);
    }
    else
    {
        allParamsVec_.resize(nParams_);

        like_.setModelCosmoParams(cosmoParams_);

        PlanckLikeFastFunc *f = new PlanckLikeFastFunc(&like_);
        func_ = f;

        PlanckLikeFastErrorFunc *errorFunc = new PlanckLikeFastErrorFunc;
        errorFunc_ = errorFunc;

        std::stringstream fileName;
        fileName << "planck_fast_" << cosmoParams_->name() << ".dat";

        layg_ = new LearnAsYouGo(cosmoParamsVec_.size() + 1, 1, *f, *errorFunc, minCount, precision, fileName.str().c_str());

        layg_->logIntoFile("planck_like_fast_log");
    }
}

PlanckLikeFast::~PlanckLikeFast()
{
    delete layg_;
    if(highT_ && !highLikeLite_)
    {
        delete (PlanckLikeFastClFunc*) func_;
        delete (PlanckLikeFastClErrorFunc*) errorFunc_;
    }
    else
    {
        delete (PlanckLikeFastFunc*) func_;
        delete (PlanckLikeFastErrorFunc*) errorFunc_;
    }

    if(logError_) outError_.close();
}

double
PlanckLikeFast::doCalculation(double* params, int nPar, bool exact)
{
    check(nPar == nParams_, "");

    const int nCosmoParams = cosmoParamsVec_.size();
    check(nCosmoParams > 0, "");
    check(nCosmoParams <= nParams_, "");

    for(int i = 0; i < nCosmoParams; ++i)
        cosmoParamsVec_[i] = params[i];

    const bool success = cosmoParams_->setAllParameters(cosmoParamsVec_);

    double error1Sigma = 0, error2Sigma = 0, errorMean = 0, errorVar = 0;

    if(!success)
    {
        const double res = 1e20 * (1.0 + rand_.generate());

        if(logError_)
        {
            for(int i = 0; i < nPar; ++i)
                outError_ << params[i] << '\t';
            outError_ << res << '\t' << error1Sigma << '\t' << error2Sigma << '\t' << errorMean << '\t' << errorVar << std::endl;
        }

        return res;
    }

    double res;

    if(highT_ && !highLikeLite_)
    {
        if(exact)
            layg_->evaluateExact(cosmoParamsVec_, &res_);
        else
            layg_->evaluate(cosmoParamsVec_, &res_, &error1Sigma, &error2Sigma, &errorMean, &errorVar);

        const int lSize = lList_.size();

        const bool usePol = (lowP_ || highP_ || lensingP_);
        const bool useBB = lowP_;
        const bool useLensing = lensingT_;

        int nCl = 1;
        if(usePol) nCl += 2;
        if(useBB) ++nCl;
        if(useLensing) ++nCl;

        check(res_.size() == nCl * lSize, "");

        std::vector<double>::const_iterator it = res_.begin();

        check(tt_.size() == lSize, "");
        for(int i = 0; i < lSize; ++i)
        {
            check(it < res_.end(), "");
            const double l = lList_[i];
            (*const_cast<std::vector<double>* >(&tt_))[i] = *(it++) * l * (l + 1);
        }

        if(usePol)
        {
            check(ee_.size() == lSize, "");
            check(te_.size() == lSize, "");

            for(int i = 0; i < lSize; ++i)
            {
                check(it < res_.end(), "");
                const double l = lList_[i];
                (*const_cast<std::vector<double>* >(&ee_))[i] = *(it++) * l * (l + 1);
            }

            for(int i = 0; i < lSize; ++i)
            {
                check(it < res_.end(), "");
                const double l = lList_[i];
                (*const_cast<std::vector<double>* >(&te_))[i] = *(it++) * l * (l + 1);
            }
        }

        if(useBB)
        {
            check(bb_.size() == lSize, "");

            for(int i = 0; i < lSize; ++i)
            {
                check(it < res_.end(), "");
                const double l = lList_[i];
                (*const_cast<std::vector<double>* >(&bb_))[i] = *(it++) * l * (l + 1);
            }
        }

        if(useLensing)
        {
            check(pp_.size() == lSize, "");

            for(int i = 0; i < lSize; ++i)
            {
                check(it < res_.end(), "");
                const double l = lList_[i];
                (*const_cast<std::vector<double>* >(&pp_))[i] = *(it++) * l * l * l * l;
            }
        }

        check(it == res_.end(), "");

        Math::CubicSpline cs(lList_, tt_);
        const int lMax = like_.getLMax();
        check(clTT_.size() == lMax + 1, "");
        for(int l = 2; l <= lMax; ++l)
            (*const_cast<std::vector<double>* >(&clTT_))[l] = cs.evaluate(double(l)) / (l * (l + 1));

        if(usePol)
        {
            Math::CubicSpline cs1(lList_, ee_);
            check(clEE_.size() == lMax + 1, "");
            for(int l = 2; l <= lMax; ++l)
                (*const_cast<std::vector<double>* >(&clEE_))[l] = cs1.evaluate(double(l)) / (l * (l + 1));

            Math::CubicSpline cs2(lList_, ee_);
            check(clTE_.size() == lMax + 1, "");
            for(int l = 2; l <= lMax; ++l)
                (*const_cast<std::vector<double>* >(&clTE_))[l] = cs2.evaluate(double(l)) / (l * (l + 1));
        }

        if(useBB)
        {
            Math::CubicSpline cs1(lList_, bb_);
            check(clBB_.size() == lMax + 1, "");
            for(int l = 2; l <= lMax; ++l)
                (*const_cast<std::vector<double>* >(&clBB_))[l] = cs1.evaluate(double(l)) / (l * (l + 1));
        }

        if(useLensing)
        {
            Math::CubicSpline cs1(lList_, pp_);
            check(clPP_.size() == lMax + 1, "");
            for(int l = 2; l <= lMax; ++l)
                (*const_cast<std::vector<double>* >(&clPP_))[l] = cs1.evaluate(double(l)) / (l * l * l * l);
        }

        like_.setCls(&clTT_, (usePol ? &clEE_ : NULL), (usePol ? &clTE_ : NULL), (useBB ? &clBB_ : NULL), (useLensing ? &clPP_ : NULL)); 

        int currentParams = nCosmoParams;

        like_.setAPlanck(params[currentParams++]);

        for(int i = 0; i < extraParams_.size(); ++i)
            extraParams_[i] = params[currentParams++];

        like_.setHighExtraParams(extraParams_);

        if(highP_)
            like_.setAPol(params[currentParams++]);

        check(currentParams == nParams_, "");

        res = like_.likelihood();
    }
    else
    {
        check(nParams_ == allParamsVec_.size(), "");

        for(int i = 0; i < nParams_; ++i)
            allParamsVec_[i] = params[i];

        if(exact)
            layg_->evaluateExact(allParamsVec_, &res_);
        else
            layg_->evaluate(allParamsVec_, &res_, &error1Sigma, &error2Sigma, &errorMean, &errorVar);

        check(res_.size() == 1, "");

        res = res_[0];
    }

    if(logError_)
    {
        for(int i = 0; i < nPar; ++i)
            outError_ << params[i] << '\t';
        outError_ << res << '\t' << error1Sigma << '\t' << error2Sigma << '\t' << errorMean << '\t' << errorVar << std::endl;
    }

    return res;
}

#else

PlanckLikeFast::PlanckLikeFast(CosmologicalParams* params, bool useCommander, bool useCamspec, bool useLensing, bool usePolarization, bool useActSpt, bool includeTensors, double kPerDecade, double precision, unsigned long minCount) : cosmoParams_(params), useCommander_(useCommander), useCamspec_(useCamspec), useLensing_(useLensing), usePol_(usePolarization), useActSpt_(useActSpt), like_(useCommander, useCamspec, useLensing, usePolarization, useActSpt, includeTensors, kPerDecade, false), logError_(false), rand_(std::time(0), 0, 1)
{
    check(useCommander_ || useCamspec_ || useLensing_ || usePol_ || useActSpt_, "at least one likelihood must be used");

    cosmoParams_->getAllParameters(cosmoParamsVec_);
    check(!cosmoParamsVec_.empty(), "");

    nParams_ = cosmoParamsVec_.size();
    if(useCamspec_)
        nParams_ += 14;
    if(useActSpt_)
        nParams_ += 24;

    lMax_ = like_.getLMax();
    cmb_.preInitialize(lMax_ + 1000, false, true, includeTensors, lMax_ + 1000, kPerDecade);
    cmb_.initialize(*cosmoParams_, true, usePol_, true);
    cmb_.getLList(lList_);

    std::vector<double>::iterator it = std::lower_bound(lList_.begin(), lList_.end(), double(lMax_));
    check(it != lList_.end(), "the l list is too short");
    ++it;
    if(it != lList_.end())
        lList_.erase(it, lList_.end());

    PlanckLikeFastFunc* f = new PlanckLikeFastFunc(cosmoParams_, &like_, &cmb_, lList_, useCommander_, useLensing_, usePol_);
    func_ = f;

    PlanckLikeFastErrorFunc* errorFunc = new PlanckLikeFastErrorFunc(&like_, lList_, useCamspec_, useActSpt_);
    errorFunc_ = errorFunc;

    std::stringstream fileName;
    fileName << "planck_fast_" << cosmoParams_->name() << ".dat";

    layg_ = new LearnAsYouGo(cosmoParamsVec_.size(), lList_.size() + 3, *f, *errorFunc, minCount, precision, fileName.str().c_str());

    layg_->logIntoFile("planck_like_fast_log");

    cl_.resize(lList_.size(), 0);
    clTT_.resize(lMax_ + 1, 0);
}

PlanckLikeFast::~PlanckLikeFast()
{
    delete layg_;
    delete (PlanckLikeFastFunc*) func_;
    delete (PlanckLikeFastErrorFunc*) errorFunc_;

    if(logError_) outError_.close();
}

double
PlanckLikeFast::doCalculation(double* params, int nPar, bool exact)
{
    check(nPar == nParams_, "");

    const int nCosmoParams = cosmoParamsVec_.size();
    check(nCosmoParams > 0, "");
    check(nCosmoParams <= nParams_, "");

    for(int i = 0; i < nCosmoParams; ++i)
        cosmoParamsVec_[i] = params[i];

    const bool success = cosmoParams_->setAllParameters(cosmoParamsVec_);

    double error1Sigma = 0, error2Sigma = 0, errorMean = 0, errorVar = 0;

    if(!success)
    {
        const double res = 1e20 * (1.0 + rand_.generate());

        if(logError_)
        {
            for(int i = 0; i < nPar; ++i)
                outError_ << params[i] << '\t';
            outError_ << res << '\t' << error1Sigma << '\t' << error2Sigma << '\t' << errorMean << '\t' << errorVar << std::endl;
        }

        return res;
    }

    if(exact)
        layg_->evaluateExact(cosmoParamsVec_, &res_);
    else
        layg_->evaluate(cosmoParamsVec_, &res_, &error1Sigma, &error2Sigma, &errorMean, &errorVar);

    const int lSize = lList_.size();
    double res = res_[lSize] + res_[lSize + 1] + res_[lSize + 2];

    if(!useCamspec_ && !useActSpt_)
        return res;

    int currentParams = nCosmoParams;

    check(cl_.size() == lList_.size(), "");
    for(int i = 0; i < cl_.size(); ++i)
    {
        const double l = lList_[i];
        cl_[i] = res_[i] * l * (l + 1);
    }

    Math::CubicSpline cs(lList_, cl_);

    check(clTT_.size() == lMax_ + 1, "");
    for(int l = 2; l <= lMax_; ++l)
        clTT_[l] = cs.evaluate(double(l)) / (l * (l + 1));

    like_.setCls(&clTT_);

    if(useCamspec_)
    {
        like_.setCamspecExtraParams(params[currentParams], params[currentParams + 1], params[currentParams + 2], params[currentParams + 3], params[currentParams + 4], params[currentParams + 5], params[currentParams + 6], params[currentParams + 7], params[currentParams + 8], params[currentParams + 9], params[currentParams + 10], params[currentParams + 11], params[currentParams + 12], params[currentParams + 13]);
        currentParams += 14;

        res += like_.camspecLike();
    }

    if(useActSpt_)
    {
        like_.setActSptExtraParams(params[currentParams], params[currentParams + 1], params[currentParams + 2], params[currentParams + 3], params[currentParams + 4], params[currentParams + 5], params[currentParams + 6], params[currentParams + 7], params[currentParams + 8], params[currentParams + 9], params[currentParams + 10], params[currentParams + 11], params[currentParams + 12], params[currentParams + 13], params[currentParams + 14], params[currentParams + 15], params[currentParams + 16], params[currentParams + 17], params[currentParams + 18], params[currentParams + 19], params[currentParams + 20], params[currentParams + 21], params[currentParams + 22], params[currentParams + 23]);
        currentParams += 24;

        res += like_.actSptLike();
    }

    check(currentParams == nParams_, "");

    if(logError_)
    {
        for(int i = 0; i < nPar; ++i)
            outError_ << params[i] << '\t';
        outError_ << res << '\t' << error1Sigma << '\t' << error2Sigma << '\t' << errorMean << '\t' << errorVar << std::endl;
    }

    return res;
}

#endif

void
PlanckLikeFast::setPrecision(double p)
{
    check(p > 0, "");
    layg_->setPrecision(p);
}
 
void
PlanckLikeFast::logError(const char* fileNameBase)
{
    check(!logError_, "this function has already been called");
    std::stringstream str;
    str << fileNameBase;
    if(CosmoMPI::create().numProcesses() > 1)
        str << "_" << CosmoMPI::create().processId();
    str << ".txt";

    outError_.open(str.str().c_str());
    if(!outError_)
    {
        StandardException exc;
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into file " << str.str();
        exc.set(exceptionStr.str());
        throw exc;
    }

    outError_ << std::setprecision(10);
    logError_ = true;
}

