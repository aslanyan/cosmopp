#include <fstream>

#include <macros.hpp>
#include <cmb.hpp>
#include <matter_likelihood.hpp>
#include <test_matter_likelihood.hpp>

std::string
TestMatterLikelihood::name() const
{
    return std::string("MATTER LIKELIHOOD TESTER");
}

unsigned int
TestMatterLikelihood::numberOfSubtests() const
{
    return 7;
}

void
TestMatterLikelihood::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < numberOfSubtests(), "invalid index " << i);

    const double h = 0.6704;
    const double omBH2 = 0.022032;
    const double omCH2 = 0.12038;
    const double tau = 0.0925;
    const double ns = 0.9619;
    const double as = 2.2154e-9;
    const double pivot = 0.05;

    const double r = 0.2;
    const double nt = 0;

    const double nEff = 3.046; 
    const int nMassive = 1;
    const double sumMNu = 0.0;

    LambdaCDMParams paramsLCDM(omBH2, omCH2, h, tau, ns, as, pivot);

    LambdaCDMParams paramsFid(0.0224, 0.274 * 0.7 * 0.7 - 0.0224, 0.7, 0.085, 0.95, 2.2e-9, pivot);

    CMB cmb;

    int lMax = 3000;
    cmb.preInitialize(lMax);
    cmb.initialize(paramsLCDM, true, false, false, true, 0.57);

    Math::TableFunction<double, double> mpk, mpkFid;
    cmb.getMatterPs(0.57, &mpk);

    cmb.initialize(paramsFid, true, false, false, true, 0.57);
    cmb.getMatterPs(0.57, &mpkFid);

    double khMin = 0, khMax = 100;
    if(i >= 4)
        khMax = 0.7;

    if(i > 1 && i < 4 || i > 4)
    {
        khMin = 0.03;
        khMax = 0.12;
    }

    MatterLikelihood ml((i < 4 ? "data/boss_dr11_pk.dat" : "data/boss_dr9_pk.dat"), (i < 4 ? "data/boss_dr11_pk_cov.dat" : "data/boss_dr9_pk_cov.dat"), i < 4, khMin, khMax);
    
    if(i > 0)
        ml.useScaling(paramsFid, 0.57);

    res = ((i == 3 || i == 6) ? ml.calculateLin(mpk, paramsLCDM) : ml.calculate(mpk, paramsLCDM));
    expected = 173.199;
    subTestName = "unscaled";

    if(i == 1)
    {
        expected = 191.083;
        subTestName = "scaled";
    }

    if(i == 2)
    {
        expected = 19.0294;
        subTestName = "scaled_range";
    }

    if(i == 3)
    {
        expected = 13.2616;
        subTestName = "scaled_range_lin";
    }

    if(i == 4)
    {
        expected = 1041.76;
        subTestName = "dr9_scaled";
    }

    if(i == 5)
    {
        expected = 35.7087;
        subTestName = "dr9_scaled_range";
    }

    if(i == 6)
    {
        expected = 28.4211;
        subTestName = "dr9_scaled_range_lin";
    }

    /*
    if(i == 1)
    {
        std::ofstream out("test_files/matter_pk_b2.txt");

        const double b2Min = 0, b2Max = 8;
        const int n = 100;
        const double b2Delta = (b2Max - b2Min) / n;
        for(int j = 0; j <= n; ++j)
        {
            const double b2 = b2Min + j * b2Delta;
            const double pk = ml.calculateLin(mpkFid, paramsFid, b2);
            out << b2 << '\t' << pk << std::endl;
        }
        out.close();
    }
    */
}
