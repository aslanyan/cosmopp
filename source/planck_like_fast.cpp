#include <string>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <cubic_spline.hpp>
#include <planck_like_fast.hpp>

namespace
{

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
        cosmoParams_->setAllParameters(x);
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
            like_->setCamspecExtraParams(100, 50, 100, 10, 30, 5, 0.9, 0.4, 0.5, 1, 1, 0.5, 5, 0.5);
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

} // namespace

PlanckLikeFast::PlanckLikeFast(CosmologicalParams* params, bool useCommander, bool useCamspec, bool useLensing, bool usePolarization, bool useActSpt, bool includeTensors, double kPerDecade, double precision, unsigned long minCount) : cosmoParams_(params), useCommander_(useCommander), useCamspec_(useCamspec), useLensing_(useLensing), usePol_(usePolarization), useActSpt_(useActSpt), like_(useCommander, useCamspec, useLensing, usePolarization, useActSpt, includeTensors, kPerDecade, false)
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
    cmb_.initialize(*cosmoParams_, true, usePol_, useLensing_);
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

    if(exact)
        layg_->evaluateExact(cosmoParamsVec_, &res_);
    else
        layg_->evaluate(cosmoParamsVec_, &res_);

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

    return res;
}

void
PlanckLikeFast::setPrecision(double p)
{
    check(p > 0, "");
    layg_->setPrecision(p);
}

