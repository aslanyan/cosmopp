#include <cosmo_mpi.hpp>

#include <fstream>
#include <sstream>
#include <memory>
#include <iomanip>

#include <macros.hpp>
#include <planck_like.hpp>
#include <mn_scanner.hpp>
#include <markov_chain.hpp>
#include <numerics.hpp>
#include <modecode.hpp>
#include <progress_meter.hpp>

namespace
{

class ModeCodeParamsUCMH : public LambdaCDMParams
{
public:
    ModeCodeParamsUCMH(double omBH2, double omCH2, double h, double tau, double kPivot, double NPivot, int potentialChoice, bool slowRollEnd, bool eternalInflOK, double kMin = 8e-7, double kMax = 1.2, int nPoints = 500) : LambdaCDMParams(omBH2, omCH2, h, tau, 1.0, 1.0, kPivot)
    {
        ModeCode::initialize(potentialChoice, kPivot, NPivot, false, false, slowRollEnd, eternalInflOK, kMin, kMax, nPoints);
        vParams_.resize(ModeCode::getNumVParams());
    }

    ~ModeCodeParamsUCMH()
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
    bool setVParams(const std::vector<double>& vParams, double *badLike) { return ModeCode::calculate(vParams, badLike); }

    virtual const Math::RealFunction& powerSpectrum() const { return ModeCode::getScalarPs(); }
    virtual const Math::RealFunction& powerSpectrumTensor() const { return ModeCode::getTensorPs(); }

    virtual void getAllParameters(std::vector<double>& v) const
    {
        check(vParams_.size() == ModeCode::getNumVParams(), "");
        v.resize(4 + vParams_.size());
        v[0] = omBH2_;
        v[1] = omCH2_;
        v[2] = h_;
        v[3] = tau_;
        for(int i = 0; i < vParams_.size(); ++i)
            v[4 + i] = vParams_[i];
    }

    virtual bool setAllParameters(const std::vector<double>& v, double *badLike)
    {
        check(v.size() == 4 + ModeCode::getNumVParams(), "");

        output_screen1("Param values:");
        for(int i = 0; i < v.size(); ++i)
            output_screen_clean1(std::setprecision(20) << "\t" << v[i]);
        output_screen_clean1(std::endl);

        setBaseParams(v[0], v[1], v[2], v[3]);

        check(vParams_.size() == ModeCode::getNumVParams(), "");
        for(int i = 0; i < vParams_.size(); ++i)
            vParams_[i] = v[4 + i];

        const bool res = setVParams(vParams_, badLike);
        output_screen1("N_piv = " << ModeCode::getNPivot() << std::endl);
        output_screen1("Result = " << res << std::endl);
        return res;
    }

private:
    std::vector<double> vParams_;
};

} // namespace

int main(int argc, char *argv[])
{
    try {
        bool ucmhLim = false;
        if(argc > 1 && std::string(argv[1]) == std::string("ucmh"))
            ucmhLim = true;

        std::string root = "slow_test_files/mn_ucmh";

#ifdef COSMO_PLANCK_15
        PlanckLikelihood like(true, true, true, true, true, false, false, true, 500);
        MnScanner mn(10, like, 500, root);
        const int nPar = 10;
#else
        PlanckLikelihood like(true, true, false, true, false, true, 500);
        MnScanner mn(23, like, 500, root);
        const int nPar = 23;
#endif

        //model 1
        //const bool slowRollEnd = true;
        //const bool eternalInflOK = false;

        //model 2
        const bool slowRollEnd = false;
        const bool eternalInflOK = true;
        ModeCodeParamsUCMH modelParams(0.02, 0.1, 0.7, 0.1, 0.002, 55, 12, slowRollEnd, eternalInflOK, 5e-6, 1.2, 500);

        if(ucmhLim)
        {
            output_screen("Adding UCMH limits!" << std::endl);
            ModeCode::addKValue(10, 0, 1e-6, 0, 1e10);
            ModeCode::addKValue(1e3, 0, 1e-7, 0, 1e10);
            ModeCode::addKValue(1e6, 0, 1e-7, 0, 1e10);
            ModeCode::addKValue(1e9, 0, 1e-2, 0, 1e10);
        }
        else
        {
            output_screen("No UCMH limits! To add these limits specify \"ucmh\" as the first argument." << std::endl);
        }

        like.setModelCosmoParams(&modelParams);

        mn.setParam(0, "ombh2", 0.02, 0.025);
        mn.setParam(1, "omch2", 0.1, 0.2);
        mn.setParam(2, "h", 0.55, 0.85);
        mn.setParam(3, "tau", 0.02, 0.20);
        //mn.setParam(4, "v_1", -10, -1);
        mn.setParam(4, "v_1", 0, 0.1);
        mn.setParam(5, "v_2", -0.1, 0.1);
        mn.setParam(6, "v_3", -0.1, 0.1);
        mn.setParam(7, "v_4", -0.1, 0.1);
        mn.setParam(8, "v_5", -12, -7);

#ifdef COSMO_PLANCK_15
        mn.setParamGauss(9, "A_planck", 1.0, 0.0025);
#else
        mn.setParam(9, "A_ps_100", 0, 360);
        mn.setParam(10, "A_ps_143", 0, 270);
        mn.setParam(11, "A_ps_217", 0, 450);
        mn.setParam(12, "A_cib_143", 0, 20);
        mn.setParam(13, "A_cib_217", 0, 80);
        mn.setParam(14, "A_sz", 0, 10);
        mn.setParam(15, "r_ps", 0.0, 1.0);
        mn.setParam(16, "r_cib", 0.0, 1.0);
        mn.setParam(17, "n_Dl_cib", -2, 2);
        mn.setParam(18, "cal_100", 0.98, 1.02);
        mn.setParam(19, "cal_127", 0.95, 1.05);
        mn.setParam(20, "xi_sz_cib", 0, 1);
        mn.setParam(21, "A_ksz", 0, 10);
        mn.setParam(22, "Bm_1_1", -20, 20);
#endif

        mn.run(true);
        
        if(!CosmoMPI::create().isMaster())
            return 0;

        MarkovChain chain("slow_test_files/mn_ucmh.txt");

        std::vector<MarkovChain::Element*> container;
        chain.getRange(container, 1.0, 0.0);

        std::ofstream outParamLimits("slow_test_files/mn_ucmh_param_limits.txt");
        for(int i = 0; i < nPar; ++i)
        {
            std::string paramName = mn.getParamName(i);

            std::stringstream fileName;
            fileName << "slow_test_files/mn_ucmh_" << paramName << ".txt";
            std::auto_ptr<Posterior1D> p(chain.posterior(i));
            p->writeIntoFile(fileName.str().c_str());

            const double median = p->median();
            double lower, upper;
            p->get1SigmaTwoSided(lower, upper);
            const double sigma = (upper - lower) / 2.0;

            outParamLimits << paramName << " = " << median << "+-" << sigma << std::endl;
        }
        outParamLimits.close();
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
