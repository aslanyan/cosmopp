#include <wmap9_like.hpp>
#include <numerics.hpp>
#include <test_wmap9_like.hpp>

std::string
TestWMAP9Like::name() const
{
    return std::string("WMAP9 LIKELIHOOD TESTER");
}

unsigned int
TestWMAP9Like::numberOfSubtests() const
{
    return 1;
}

void
TestWMAP9Like::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
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

    const bool useLowlT = true, useHighlT = true, useLowlP = true, useHighlP = true, useGibbs = true, useTTBeam = true;

    WMAP9Likelihood like(useLowlT, useHighlT, useLowlP, useHighlP, useGibbs, useTTBeam);

    like.setCosmoParams(paramsLCDM);
    like.calculateCls();

    subTestName = std::string("LCDM");
    res = like.likelihood();

    check(Math::areEqual(res, like.lowlTLike() + like.highlTLike() + like.TTBeamLike() + like.lowlPLike() + like.TELike() + like.TBLike(), 1e-5), "");

    output_screen1("low-l TT chi2 = " << like.lowlTChi2() << ", det = " << like.lowlTDet() << ", total = " << like.lowlTLike() << std::endl);
    output_screen1("high-l TT like = " << like.highlTLike() << std::endl);
    output_screen1("TT beam like = " << like.TTBeamLike() << std::endl);
    output_screen1("low-l P chi2 = " << like.lowlPChi2() << ", det = " << like.lowlPDet() << ", total = " << like.lowlPLike() << std::endl);
    output_screen1("high-l TE chi2 = " << like.TEChi2() << ", det = " << like.TEDet() << ", total = " << like.TELike() << std::endl);
    output_screen1("high-l TB chi2 = " << like.TBChi2() << ", det = " << like.TBDet() << ", total = " << like.TBLike() << std::endl);

    expected = 7592.002623; // calculated by running the WMAP9 likelihood code itself on the same cl values
}
