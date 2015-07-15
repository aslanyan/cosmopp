#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#include <test_multinest_planck.hpp>
#include <mn_scanner.hpp>
#include <planck_like.hpp>
#include <markov_chain.hpp>
#include <numerics.hpp>
#include <timer.hpp>

std::string
TestMultinestPlanck::name() const
{
    return std::string("MULTINEST PLANCK LIKELIHOOD TESTER");
}

unsigned int
TestMultinestPlanck::numberOfSubtests() const
{
    return 1;
}

void
TestMultinestPlanck::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 1, "invalid index " << i);
    
    using namespace Math;

    std::string root = "slow_test_files/multinest_planck_test";
#ifdef COSMO_PLANCK_15
    PlanckLikelihood planckLike(true, true, true, false, true, false, false, false, 5);
    MnScanner mn(7, planckLike, 300, root);
#else
    PlanckLikelihood planckLike(true, true, false, true, false, false, 5);
    MnScanner mn(20, planckLike, 300, root);
#endif

    mn.setParam(0, "ombh2", 0.02, 0.025);
    mn.setParam(1, "omch2", 0.1, 0.2);
    mn.setParam(2, "h", 0.55, 0.85);
    mn.setParam(3, "tau", 0.02, 0.20);
    mn.setParam(4, "ns", 0.9, 1.1);
    mn.setParam(5, "As", 2.7, 3.5);

#ifdef COSMO_PLANCK_15
    mn.setParamGauss(6, "A_planck", 1.0, 0.0025);
#else
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
#endif

    const double pivot = 0.05;
    LambdaCDMParams par(0.022, 0.12, 0.7, 0.1, 1.0, std::exp(3.0) / 1e10, pivot);
    planckLike.setModelCosmoParams(&par);

    Timer timer("MN PLANCK");
    timer.start();
    mn.run(true);
    const unsigned long time = timer.end();
    output_screen("MN Planck took " << time / 1000000 << " seconds." << std::endl);
    
    subTestName = std::string("standard_param_limits");
    res = 1;
    expected = 1;

    if(!isMaster())
        return;

    MarkovChain chain("slow_test_files/multinest_planck_test.txt");

    const int nPoints = 1000;

#ifdef COSMO_PLANCK_15
    const double expectedMedian[6] = {0.02222, 0.1197, 0.6731, 0.078, 0.9655, 3.089};
    const double expectedSigma[6] = {0.00023, 0.0022, 0.0096, 0.019, 0.0062, 0.036};
    const int nPar = 7;
#else
    const double expectedMedian[6] = {0.02205, 0.1199, 0.673, 0.089, 0.9603, 3.089};
    const double expectedSigma[6] = {0.00028, 0.0027, 0.012, 0.013, 0.0073, 0.025};
    const int nPar = 20;
#endif

    std::ofstream outParamLimits("slow_test_files/multinest_planck_param_limits.txt");
    for(int i = 0; i < nPar; ++i)
    {
        const std::string& paramName = mn.getParamName(i);
        std::stringstream fileName;
        fileName << "slow_test_files/multinest_planck_" << paramName << ".txt";
        Posterior1D* p = chain.posterior(i);

        p->writeIntoFile(fileName.str().c_str(), nPoints);

        const double median = p->median();
        double lower, upper;
        p->get1SigmaTwoSided(lower, upper);
        const double sigma = (upper - lower) / 2.0;

        outParamLimits << paramName << " = " << median << "+-" << sigma << std::endl;
        // check the standard cosmological parameter limits
        if(i < 6)
        {
            if(std::abs(expectedMedian[i] - median) > expectedSigma[i])
            {
                output_screen("FAIL: Expected " << paramName << " median is " << expectedMedian[i] << ", the result is " << median << std::endl);
                res = 0;
            }

            if(!Math::areEqual(expectedSigma[i], sigma, 0.25))
            {
                output_screen("FAIL: Expected" << paramName << "sigma is " << expectedSigma[i] << ", the result is " << sigma << std::endl);
                res = 0;
            }
        }
    }
    outParamLimits.close();
}
