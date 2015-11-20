#include <map>
#include <memory>

#include <cosmological_params.hpp>
#include <polychord.hpp>
#include <planck_like_fast.hpp>

class LCDMRunRunParams : public LambdaCDMParams
{
public:
    LCDMRunRunParams(double omBH2, double omCH2, double h, double tau, double ns, double as, double pivot, double run = 0.0) : LambdaCDMParams(omBH2, omCH2, h, tau, ns, as, pivot, run) {}
    virtual ~LCDMRunRunParams() {}

    virtual std::string name() const { return "LCDMRunRun"; }
    virtual void getAllParameters(std::vector<double>& v) const
    {
        v.resize(8);
        v[0] = getOmBH2();
        v[1] = getOmCH2();
        v[2] = getH();
        v[3] = getTau();
        v[4] = getNs();
        v[5] = std::log(getAs() * 1e10);
        v[6] = ps_.getRun();
        v[7] = ps_.getRunRun();
    }

    virtual bool setAllParameters(const std::vector<double>& v, double *badLike)
    {
        check(v.size() >= 8, "");
        omBH2_ = v[0];
        omCH2_ = v[1];
        h_ = v[2];
        tau_ = v[3];
        ps_.setNs(v[4]);
        ps_.setAs(std::exp(v[5]) / 1e10);
        ps_.setRun(v[6]);
        ps_.setRunRun(v[7]);

        double bad = 0;

        for(std::map<double, double>::const_iterator it = kLim_.begin(); it != kLim_.end(); ++it)
        {
            if(ps_.evaluate(it->first) > it->second)
                bad += (ps_.evaluate(it->first) - it->second) / it->second;
        }

        check(bad >= 0, "");
        if(badLike)
            *badLike = bad * 1e10;

        return (bad == 0);
    }

    void addKLimit(double k, double lim) { kLim_[k] = lim; }

private:
    std::map<double, double> kLim_;
};

int main(int argc, char *argv[])
{
    try {
        bool ucmhLim = false;
        //bool useFast = false;

        for(int i = 1; i < argc; ++i)
        {
            if(std::string(argv[i]) == "ucmh")
                ucmhLim = true;
            //if(std::string(argv[i]) == "fast")
                //useFast = true;
        }

        using namespace Math;

        const double pivot = 0.05;
        LCDMRunRunParams par(0.022, 0.12, 0.7, 0.1, 1.0, std::exp(3.0) / 1e10, pivot);

        if(ucmhLim)
        {
            output_screen("Adding UCMH limits!" << std::endl);
            par.addKLimit(10, 1e-6);
            par.addKLimit(1e3, 1e-7);
            par.addKLimit(1e6, 1e-7);
            par.addKLimit(1e9, 1e-2);
        }
        else
        {
            output_screen("No UCMH limits! To add these limits specify \"ucmh\" as an argument." << std::endl);
        }

        std::stringstream root;
        root << "slow_test_files/pc_planck_run_run";
        //if(useFast)
            //root << "_fast";

        std::unique_ptr<PlanckLikelihood> planckLike;
        /*
        if(useFast)
        {
            output_screen("Using the fast version of Planck likelihood!" << std::endl);
            PlanckLikeFast *like = new PlanckLikeFast(&par, true, true, false, true, false, false, 5, 0.4, 50000);
            std::string errorLogRoot = "slow_test_files/pc_planck_run_run_fast_error_log";
            like->logError(errorLogRoot.c_str());
            planckLike.reset(like);
        }
        else
        */
        //{
            //output_screen("Using the regular version of Planck likelihood! To use the fast version specify \"fast\" as an argument." << std::endl);
#ifdef COSMO_PLANCK_15
            PlanckLikelihood *like = new PlanckLikelihood(true, true, true, false, true, false, false, false, 100, true);
            like->setModelCosmoParams(&par);
            planckLike.reset(like);
            PolyChord pc(9, *planckLike, 300, root.str(), 4);
#else
            PlanckLikelihood *like = new PlanckLikelihood(true, true, false, true, false, false, 100);
            like->setModelCosmoParams(&par);
            planckLike.reset(like);
            PolyChord pc(22, *planckLike, 300, root.str(), 4);
#endif
        //}


        pc.setParam(0, "ombh2", 0.02, 0.025, 1);
        pc.setParam(1, "omch2", 0.1, 0.2, 1);
        pc.setParam(2, "h", 0.55, 0.85, 1);
        pc.setParam(3, "tau", 0.02, 0.20, 1);
        pc.setParam(4, "ns", 0.9, 1.1, 2);
        pc.setParam(5, "As", 2.7, 3.5, 2);
        pc.setParam(6, "nrun", -1.0, 1.0, 2);
        pc.setParam(7, "nrunrun", -1.0, 1.0, 2);

#ifdef COSMO_PLANCK_15
        pc.setParamGauss(8, "A_planck", 1.0, 0.0025, 3);
#else
        pc.setParam(8, "A_ps_100", 0, 360, 3);
        pc.setParam(9, "A_ps_143", 0, 270, 3);
        pc.setParam(10, "A_ps_217", 0, 450, 3);
        pc.setParam(11, "A_cib_143", 0, 20, 3);
        pc.setParam(12, "A_cib_217", 0, 80, 3);
        pc.setParam(13, "A_sz", 0, 10, 3);
        pc.setParam(14, "r_ps", 0.0, 1.0, 3);
        pc.setParam(15, "r_cib", 0.0, 1.0, 3);
        pc.setParam(16, "n_Dl_cib", -2, 2, 3);
        pc.setParam(17, "cal_100", 0.98, 1.02, 3);
        pc.setParam(18, "cal_127", 0.95, 1.05, 3);
        pc.setParam(19, "xi_sz_cib", 0, 1, 3);
        pc.setParam(20, "A_ksz", 0, 10, 3);
        pc.setParam(21, "Bm_1_1", -20, 20, 3);
#endif

        if(ucmhLim)
        {
            const std::vector<double> fracs{0.4, 0.5, 0.1};
            pc.setParameterHierarchy(fracs);
        }
        else
        {
            const std::vector<double> fracs{0.80, 0.18, 0.02};
            pc.setParameterHierarchy(fracs);
        }

        Timer timer("MN PLANCK RUN_RUN");
        timer.start();
        pc.run(true);
        const unsigned long time = timer.end();
        output_screen("MN Planck run_run took " << time / 1000000 << " seconds." << std::endl);
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
