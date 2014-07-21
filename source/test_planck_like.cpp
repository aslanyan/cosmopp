#include <planck_like.hpp>
#include <test_planck_like.hpp>

std::string
TestPlanckLike::name() const
{
    return std::string("PLANCK LIKELIHOOD TESTER");
}

unsigned int
TestPlanckLike::numberOfSubtests() const
{
    return 2;
}

void
TestPlanckLike::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
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
    LCDMWithTensorParams paramsLCDMTens(omBH2, omCH2, h, tau, ns, as, pivot, r, nt, pivot); 

    PlanckLikelihood like(true, true, true, true, false, bool(i));

    if(i == 0)
        like.setCosmoParams(paramsLCDM);
    else
        like.setCosmoParams(paramsLCDMTens);
    like.setCamspecExtraParams(153, 54.9, 55.8, 4, 55.5, 4, 0.91, 0.63, 0.6, 1, 1, 0.1, 1, 0.3);

    like.calculateCls();

    subTestName = std::string("LCDM");
    if(i == 1)
        subTestName = std::string("LCDM+Tensor");
    res = like.likelihood();
    expected = (i == 0 ? 10038.47 : 10048.5); // calculated by running the Planck likelihood code itself on the same cl values
}

