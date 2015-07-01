#include <cosmo_mpi.hpp>

#include <fstream>
#include <sstream>
#include <string>
#include <memory>

#include <macros.hpp>
#include <planck_like.hpp>
#include <mn_scanner.hpp>
#include <markov_chain.hpp>
#include <numerics.hpp>
#include <modecode.hpp>
#include <progress_meter.hpp>

int main(int argc, char *argv[])
{
    try {
        bool ucmhLim = false;
        if(argc > 1 && std::string(argv[1]) == std::string("ucmh"))
            ucmhLim = true;

        PlanckLikelihood like(true, true, false, true, false, true, 500);
        ModeCodeCosmologicalParams modelParams(0.02, 0.1, 0.7, 0.1, 0.002, 55, 12, false, true, false, false, 8e-7, 1.2, 500);

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

        std::string root = "slow_test_files/mn_ucmh";
        MnScanner mn(24, like, 500, root);

        mn.setParam(0, "ombh2", 0.02, 0.025);
        mn.setParam(1, "omch2", 0.1, 0.2);
        mn.setParam(2, "h", 0.55, 0.85);
        mn.setParam(3, "tau", 0.02, 0.20);
        mn.setParam(4, "NPivot", 35, 90);
        //mn.setParam(4, "NPivot", 20, 90);
        mn.setParam(5, "v_1", -10, -1);
        mn.setParam(6, "v_2", -0.1, 0.1);
        mn.setParam(7, "v_3", -0.1, 0.1);
        //mn.setParam(8, "v_4", 0.0, 0.0);
        mn.setParam(8, "v_4", -0.1, 0.1);
        mn.setParam(9, "v_5", -12, -9);

        mn.setParam(10, "A_ps_100", 0, 360);
        mn.setParam(11, "A_ps_143", 0, 270);
        mn.setParam(12, "A_ps_217", 0, 450);
        mn.setParam(13, "A_cib_143", 0, 20);
        mn.setParam(14, "A_cib_217", 0, 80);
        mn.setParam(15, "A_sz", 0, 10);
        mn.setParam(16, "r_ps", 0.0, 1.0);
        mn.setParam(17, "r_cib", 0.0, 1.0);
        mn.setParam(18, "n_Dl_cib", -2, 2);
        mn.setParam(19, "cal_100", 0.98, 1.02);
        mn.setParam(20, "cal_127", 0.95, 1.05);
        mn.setParam(21, "xi_sz_cib", 0, 1);
        mn.setParam(22, "A_ksz", 0, 10);
        mn.setParam(23, "Bm_1_1", -20, 20);

        mn.run(true);
        
        if(!CosmoMPI::create().isMaster())
            return 0;

        MarkovChain chain("slow_test_files/mn_ucmh.txt");

        std::vector<MarkovChain::Element*> container;
        chain.getRange(container, 1.0, 0.0);

        std::ofstream outParamLimits("slow_test_files/mn_ucmh_param_limits.txt");
        for(int i = 0; i < 24; ++i)
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
