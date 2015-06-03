#include <fstream>
#include <sstream>

#include <macros.hpp>
#include <planck_like.hpp>
#include <mn_scanner.hpp>
#include <markov_chain.hpp>
#include <numerics.hpp>
#include <modecode.hpp>
#include <progress_meter.hpp>
#include <test_mn_modecode.hpp>

std::string
TestMNModeCode::name() const
{
    return std::string("MULTINEST MODECODE TESTER");
}

unsigned int
TestMNModeCode::numberOfSubtests() const
{
    return 1;
}

/*
class MSqPhiSqCosmologicalParams : public LambdaCDMParams
{
public:
    MSqPhiSqCosmologicalParams(double omBH2, double omCH2, double h, double tau, double kPivot, double NPivot, double logMSq, bool instReheat = true) : LambdaCDMParams(omBH2, omCH2, h, tau, 1.0, 1.0, kPivot), vParams_(1, logMSq)
    {
        ModeCode::initialize(1, kPivot, NPivot, instReheat, 1e-6, 1.0, 500);
        ModeCode::calculate(vParams_);
    }

    ~MSqPhiSqCosmologicalParams()
    {
    }

    void setBaseParams(double omBH2, double omCH2, double h, double tau)
    {
        omBH2_ = omBH2;
        omCH2_ = omCH2;
        h_ = h;
        tau_ = tau;
    }

    void setNPivot(double NPivot) { ModeCode::setNPivot(NPivot); }
    void setLogMSq(double logMSq) { check(vParams_.size() == 1, ""); vParams_[0] = logMSq; ModeCode::calculate(vParams_); }

    virtual const Math::RealFunction& powerSpectrum() const { return ModeCode::getScalarPs(); }
    virtual const Math::RealFunction& powerSpectrumTensor() const { return ModeCode::getTensorPs(); }

private:
    std::vector<double> vParams_;
};
*/

/*
class TestMNModeCodeLike : public Math::LikelihoodFunction
{
public:
    enum PSTYPE { LCDM = 0, MSQPHISQ, PSTYPE_MAX };

public:
    TestMNModeCodeLike(PSTYPE psType = LCDM) : psType_(psType), kPivot_(0.05), like_(true, true, false, true, false, false), vParams_(1, 0)
    {
        check(psType_ >= 0 && psType_ < PSTYPE_MAX, "invalid ps type " << psType_);

        if(psType_ == MSQPHISQ)
            params_ = new ModeCodeCosmologicalParams(0.02, 0.1, 0.7, 0.1, kPivot_, 55, 1, true, true);
        else
            params_ = NULL;
    }
    ~TestMNModeCodeLike()
    {
        if(params_)
            delete params_;
    }

    virtual double calculate(double* params, int nParams)
    {
        check(nParams == 20, "");
        if(psType_ == LCDM)
        {
            check(!params_, "");
            params_ = new LambdaCDMParams(params[0], params[1], params[2], params[3], params[4], std::exp(params[5]) / 1e10, kPivot_);
        }
        else
        {
            check(psType_ == MSQPHISQ, "");
            check(params_, "");
            ModeCodeCosmologicalParams* p = (ModeCodeCosmologicalParams*) params_;
            p->setBaseParams(params[0], params[1], params[2], params[3]);
            p->setNPivot(params[4]);
            
            check(vParams_.size() == 1, "");
            vParams_[0] = params[5];
            p->setVParams(vParams_);
        }

        like_.setCosmoParams(*params_);
        like_.setCamspecExtraParams(params[6], params[7], params[8], params[9], params[10], params[11], params[12], params[13], params[14], params[15], params[16], params[17], params[18], params[19]);

        const double like = like_.likelihood();

        if(psType_ == LCDM)
        {
            delete params_;
            params_ = NULL;
        }

        return like;
    }

private:
    const double kPivot_;
    PSTYPE psType_;
    CosmologicalParams* params_;
    PlanckLikelihood like_;
    std::vector<double> vParams_;
};
*/

void
TestMNModeCode::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 1, "invalid index " << i);

    PlanckLikelihood like(true, true, false, true, false, false);
    ModeCodeCosmologicalParams modelParams(0.02, 0.1, 0.7, 0.1, 0.05, 55, 1, true, false, true, 8e-7, 1.2, 500);
    like.setModelCosmoParams(&modelParams);

    std::string root = "slow_test_files/mn_modecode_test";
    MnScanner mn(20, like, 300, root);

    mn.setParam(0, "ombh2", 0.02, 0.025);
    mn.setParam(1, "omch2", 0.1, 0.2);
    mn.setParam(2, "h", 0.55, 0.85);
    mn.setParam(3, "tau", 0.02, 0.20);
    mn.setParam(4, "NPivot", 20, 90);
    mn.setParam(5, "log_msq", -13.5, -8);

    mn.setParam(6, "A_ps_100", 0, 360);
    mn.setParam(7, "A_ps_143", 0, 270);
    mn.setParam(8, "A_ps_217", 0, 450);
    mn.setParam(9, "A_cib_143", 0, 20);
    mn.setParam(10, "A_cib_217", 0, 80);
    mn.setParam(11, "A_sz", 0, 10);
    mn.setParam(12, "r_ps", 0.0, 1.0);
    mn.setParam(13, "r_cib", 0.0, 1.0);
    mn.setParam(14, "n_Dl_cib", -2, 2);
    mn.setParam(15, "cal_100", 0.98, 1.02);
    mn.setParam(16, "cal_127", 0.95, 1.05);
    mn.setParam(17, "xi_sz_cib", 0, 1);
    mn.setParam(18, "A_ksz", 0, 10);
    mn.setParam(19, "Bm_1_1", -20, 20);

    mn.run(true);
    
    subTestName = std::string("msqphisq_param_limits");
    res = 1;
    expected = 1;

    if(!isMaster())
        return;

    MarkovChain chain("slow_test_files/mn_modecode_test.txt");

    //TBD
    //const double expectedMedian[6] = {0.02205, 0.1199, 0.673, 0.089, 0.9603, 3.089};
    //const double expectedSigma[6] = {0.00028, 0.0027, 0.012, 0.013, 0.0073, 0.025};
    const double expectedMedian[6] = {0.02223, 0.1183, 0.686, 0.088, 55, -10.41};
    const double expectedSigma[6] = {0.00023, 0.0015, 0.007, 0.012, 22, 0.01};

    std::vector<MarkovChain::Element*> container;
    chain.getRange(container, 1.0, 0.0);
    Posterior1D nsPost, rPost;

    const double tensorPivot = 0.002;
    const double kMin = 1e-6, kMax = 1.0;
    ModeCode::initialize(1, 0.05, 55, true, kMin, kMax, 100);
    std::vector<double> vPar(1, 0);

    output_screen("Generating the posterior distributions for ns and r..." << std::endl);
    ProgressMeter meter(container.size());
    for(unsigned long j = 0; j < container.size(); ++j)
    {
        const double NPivot = container[j]->params[4], logMsq = container[j]->params[5];
        ModeCode::setNPivot(NPivot);
        vPar[0] = logMsq;
        ModeCode::calculate(vPar);
        const double sMin = ModeCode::getScalarPs().evaluate(kMin);
        const double sMax = ModeCode::getScalarPs().evaluate(kMax);

        const double ns = 1 + (std::log(sMax) - std::log(sMin)) / (std::log(kMax) - std::log(kMin));
        const double r = ModeCode::getTensorPs().evaluate(tensorPivot) / ModeCode::getScalarPs().evaluate(tensorPivot);

        nsPost.addPoint(ns, container[j]->prob, container[j]->like);
        rPost.addPoint(r, container[j]->prob, container[j]->like);

        meter.advance();
    }

    output_screen("OK" << std::endl);

    nsPost.generate();

    output_screen("ns done!" << std::endl);
    rPost.generate();
    output_screen("r done!" << std::endl);

    std::ofstream outParamLimits("slow_test_files/mn_modecode_param_limits.txt");
    for(int i = 0; i < 22; ++i)
    {
        std::string paramName;
        
        if(i < 20)
            mn.getParamName(i);

        if(i == 20)
            paramName = "ns";
        if(i == 21)
            paramName = "r";

        std::stringstream fileName;
        fileName << "slow_test_files/mn_modecode_" << paramName << ".txt";
        Posterior1D* p = (i < 20 ? chain.posterior(i) : (i == 20 ? &nsPost : &rPost));
        p->writeIntoFile(fileName.str().c_str());

        const double median = p->median();
        double lower, upper;
        p->get1SigmaTwoSided(lower, upper);
        const double sigma = (upper - lower) / 2.0;

        outParamLimits << paramName << " = " << median << "+-" << sigma << std::endl;
        // check the standard cosmological parameter limits
        if(i < 4)
        {
            if(std::abs(expectedMedian[i] - median) > expectedSigma[i] / 2)
            {
                output_screen("FAIL: Expected " << paramName << " median is " << expectedMedian[i] << ", the result is " << median << std::endl);
                res = 0;
            }
        }

        if(i < 20)
            delete p;
    }
    outParamLimits.close();
}
