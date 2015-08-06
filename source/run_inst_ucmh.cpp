#include <cosmo_mpi.hpp>

#include <fstream>
#include <sstream>
#include <memory>
#include <iomanip>
#include <memory>

#include <macros.hpp>
#include <planck_like.hpp>
#include <mn_scanner.hpp>
#include <mcmc.hpp>
#include <markov_chain.hpp>
#include <numerics.hpp>
#include <modecode.hpp>
#include <taylor_pk.hpp>
#include <progress_meter.hpp>

namespace
{

class TaylorParamsUCMH : public LambdaCDMParams
{
public:
    TaylorParamsUCMH(double omBH2, double omCH2, double h, double tau, double kPivot, double NPivot, int potentialChoice, bool slowRollEnd, bool eternalInflOK, double kMin = 8e-7, double kMax = 1.2, int nPoints = 500, bool useClass = false) : LambdaCDMParams(omBH2, omCH2, h, tau, 1.0, 1.0, kPivot), useClass_(useClass)
    {
        if(useClass)
        {
            taylor_.reset(new TaylorPk(kPivot, kMin, kMax, 10));
            vParams_.resize(5);
        }
        else
        {
            ModeCode::initialize(potentialChoice, kPivot, NPivot, false, false, slowRollEnd, eternalInflOK, kMin, kMax, nPoints);
            vParams_.resize(ModeCode::getNumVParams());
        }
    }

    ~TaylorParamsUCMH()
    {
    }

    void addKValue(double k, double sMin = 0, double sMax = 1, double tMin = 0, double tMax = 1)
    {
        if(useClass_)
            taylor_->addKValue(k, sMin, sMax, tMin, tMax);
        else
            ModeCode::addKValue(k, sMin, sMax, tMin, tMax);
    }

    void setBaseParams(double omBH2, double omCH2, double h, double tau)
    {
        omBH2_ = omBH2;
        omCH2_ = omCH2;
        h_ = h;
        tau_ = tau;
    }

    bool setVParams(const std::vector<double>& vParams, double *badLike)
    {
        if(useClass_)
            return taylor_->calculate(vParams, badLike);

        return ModeCode::calculate(vParams, badLike);
    }

    virtual const Math::RealFunction& powerSpectrum() const
    {
        if(useClass_) return taylor_->getScalarPs();

        return ModeCode::getScalarPs();
    }

    virtual const Math::RealFunction& powerSpectrumTensor() const
    {
        if(useClass_) return taylor_->getTensorPs();

        return ModeCode::getTensorPs();
    }

    virtual void getAllParameters(std::vector<double>& v) const
    {
        check(useClass_ || vParams_.size() == ModeCode::getNumVParams(), "");
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
        check(useClass_ || v.size() == 4 + ModeCode::getNumVParams(), "");

        output_screen1("Param values:");
        //output_log("Param values:");
        for(int i = 0; i < v.size(); ++i)
        {
            output_screen_clean1(std::setprecision(20) << "\t" << v[i]);
            //output_log(std::setprecision(20) << "\t" << v[i]);
        }
        output_screen_clean1(std::endl);
        //output_log(std::endl);

        setBaseParams(v[0], v[1], v[2], v[3]);

        check(useClass_ || vParams_.size() == ModeCode::getNumVParams(), "");
        for(int i = 0; i < vParams_.size(); ++i)
            vParams_[i] = v[4 + i];

        const bool res = setVParams(vParams_, badLike);
        if(!useClass_)
        {
            output_screen1("N_piv = " << ModeCode::getNPivot() << std::endl);
        }
        output_screen1("Result = " << res << std::endl);
        return res;
    }

private:
    const bool useClass_;
    std::auto_ptr<TaylorPk> taylor_;
    std::vector<double> vParams_;
};

} // namespace

int main(int argc, char *argv[])
{
    try {
        bool ucmhLim = false;
        bool useClass = false;
        bool useMH = false;
        for(int i = 1; i < argc; ++i)
        {
            if(std::string(argv[i]) == std::string("ucmh"))
                ucmhLim = true;

            if(std::string(argv[i]) == std::string("class"))
                useClass = true;

            if(std::string(argv[i]) == std::string("mh"))
                useMH = true;
        }

        if(useClass)
        {
            output_screen("Using CLASS for calculating pk." << std::endl);
        }
        else
        {
            output_screen("Using Modecode for calculating pk. To use CLASS instead specify \"class\" as an argument." << std::endl);
        }

        if(useMH)
        {
            output_screen("Using Metropolis-Hastings sampler." << std::endl);
        }
        else
        {
            output_screen("Using MultiNest sampler. To use Metropolis-Hastings instead specify \"mh\" as an argument." << std::endl);
        }

        std::string root = (useMH ? "slow_test_files/mh_ucmh" : "slow_test_files/mn_ucmh");

        std::auto_ptr<MnScanner> mn;
        std::auto_ptr<Math::MetropolisHastings> mh;

#ifdef COSMO_PLANCK_15
        PlanckLikelihood like(true, true, true, true, true, false, false, true, 500);
        const int nPar = 10;
#else
        PlanckLikelihood like(true, true, false, true, false, true, 500);
        const int nPar = 23;
#endif

        if(useMH)
            mh.reset(new Math::MetropolisHastings(nPar, like, root));
        else
            mn.reset(new MnScanner(nPar, like, 500, root));

        const double kPivot = 0.05;

        //model 1
        //const bool slowRollEnd = true;
        //const bool eternalInflOK = false;

        //model 2
        const bool slowRollEnd = false;
        const bool eternalInflOK = true;
        TaylorParamsUCMH modelParams(0.02, 0.1, 0.7, 0.1, kPivot, 55, 12, slowRollEnd, eternalInflOK, 5e-6, 1.2, 500, useClass);

        if(ucmhLim)
        {
            output_screen("Adding UCMH limits!" << std::endl);
            modelParams.addKValue(10, 0, 1e-6, 0, 1e10);
            modelParams.addKValue(1e3, 0, 1e-7, 0, 1e10);
            modelParams.addKValue(1e6, 0, 1e-7, 0, 1e10);
            modelParams.addKValue(1e9, 0, 1e-2, 0, 1e10);
        }
        else
        {
            output_screen("No UCMH limits! To add these limits specify \"ucmh\" as an argument." << std::endl);
        }

        like.setModelCosmoParams(&modelParams);

        int nChains;
        unsigned long burnin;
        unsigned int thin;

        if(useMH)
        {
            mh->setParam(0, "ombh2", 0.02, 0.025, 0.022, 0.0003, 0.0001);
            mh->setParam(1, "omch2", 0.1, 0.2, 0.12, 0.003, 0.001);
            mh->setParam(2, "h", 0.55, 0.85, 0.68, 0.02, 0.005);
            mh->setParam(3, "tau", 0.02, 0.2, 0.1, 0.02, 0.01);
            mh->setParam(4, "v_1", 0, 0.1, 0.01, 0.005, 0.005);
            mh->setParam(5, "v_2", -0.1, 0.1, 0, 0.02, 0.02);
            mh->setParam(6, "v_3", -0.1, 0.1, 0, 0.01, 0.01);
            mh->setParam(7, "v_4", -0.1, 0.1, 0, 0.01, 0.01);
            mh->setParam(8, "v_5", -12, -8, -9, 0.5, 0.1);

#ifdef COSMO_PLANCK_15
            mh->setParamGauss(9, "A_planck", 1.0, 0.0025, 1.0, 0.002, 0.002);
#else
            mh->setParam(9, "A_ps_100", 0, 360, 100, 100, 20);
            mh->setParam(10, "A_ps_143", 0, 270, 50, 20, 2);
            mh->setParam(11, "A_ps_217", 0, 450, 100, 30, 4);
            mh->setParam(12, "A_cib_143", 0, 20, 10, 10, 1);
            mh->setParam(13, "A_cib_217", 0, 80, 30, 15, 1);
            mh->setParam(14, "A_sz", 0, 10, 5, 5, 1);
            mh->setParam(15, "r_ps", 0.0, 1.0, 0.9, 0.2, 0.02);
            mh->setParam(16, "r_cib", 0.0, 1.0, 0.4, 0.4, 0.05);
            mh->setParam(17, "n_Dl_cib", -2, 2, 0.5, 0.2, 0.02);
            mh->setParam(18, "cal_100", 0.98, 1.02, 1.0, 0.0008, 0.0001);
            mh->setParam(19, "cal_127", 0.95, 1.05, 1.0, 0.003, 0.0002);
            mh->setParam(20, "xi_sz_cib", 0, 1, 0.5, 0.6, 0.05);
            mh->setParam(21, "A_ksz", 0, 10, 5, 6, 0.5);
            mh->setParam(22, "Bm_1_1", -20, 20, 0.5, 1.0, 0.1);
#endif

            burnin = 1000;
            thin = 2;
            nChains = mh->run(100000, 1, burnin, Math::MetropolisHastings::GELMAN_RUBIN, 0.01, false);
        }
        else
        {
            mn->setParam(0, "ombh2", 0.02, 0.025);
            mn->setParam(1, "omch2", 0.1, 0.2);
            mn->setParam(2, "h", 0.55, 0.85);
            mn->setParam(3, "tau", 0.02, 0.20);
            mn->setParam(4, "v_1", 0, 0.1);
            mn->setParam(5, "v_2", -0.1, 0.1);
            mn->setParam(6, "v_3", -0.1, 0.1);
            mn->setParam(7, "v_4", -0.1, 0.1);
            mn->setParam(8, "v_5", -12, -7);

#ifdef COSMO_PLANCK_15
            mn->setParamGauss(9, "A_planck", 1.0, 0.0025);
#else
            mn->setParam(9, "A_ps_100", 0, 360);
            mn->setParam(10, "A_ps_143", 0, 270);
            mn->setParam(11, "A_ps_217", 0, 450);
            mn->setParam(12, "A_cib_143", 0, 20);
            mn->setParam(13, "A_cib_217", 0, 80);
            mn->setParam(14, "A_sz", 0, 10);
            mn->setParam(15, "r_ps", 0.0, 1.0);
            mn->setParam(16, "r_cib", 0.0, 1.0);
            mn->setParam(17, "n_Dl_cib", -2, 2);
            mn->setParam(18, "cal_100", 0.98, 1.02);
            mn->setParam(19, "cal_127", 0.95, 1.05);
            mn->setParam(20, "xi_sz_cib", 0, 1);
            mn->setParam(21, "A_ksz", 0, 10);
            mn->setParam(22, "Bm_1_1", -20, 20);
#endif

            mn->run(true);

            nChains = 1;
            burnin = 0;
            thin = 1;
        }
        
        if(!CosmoMPI::create().isMaster())
            return 0;

        MarkovChain chain(nChains, root.c_str(), burnin, thin);

        std::vector<MarkovChain::Element*> container;
        chain.getRange(container, 1.0, 0.0);

        std::ofstream outParamLimits("slow_test_files/mn_ucmh_param_limits.txt");
        for(int i = 0; i < nPar; ++i)
        {
            std::string paramName = (useMH ? mh->getParamName(i) : mn->getParamName(i));

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
