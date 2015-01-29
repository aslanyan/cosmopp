#include <cmb.hpp>
#include <test_cmb.hpp>

std::string
TestCMB::name() const
{
    return std::string("CMB TESTER");
}

unsigned int
TestCMB::numberOfSubtests() const
{
    return 1;
}

void
TestCMB::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 1, "invalid index " << i);

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
    LCDMWithTensorParams paramsLCDMTens(omBH2, omCH2, h, tau, ns, as, pivot, r, nt, pivot); 

    CMB cmb;

    int lMax = 3000;
    std::vector<double> clTT;

    switch(i)
    {
    case 0:
        cmb.preInitialize(lMax, false, true, false);
        cmb.initialize(paramsLCDM, true, false, true, false);
        cmb.getLensedCl(&clTT);

        subTestName = std::string("LCDM");
        res = clTT[10];
        expected = 47.8189;
        break;
    default:
        check(false, "");
        break;
    }
}
