#include <memory>

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
#ifdef COSMO_PLANCK_15
    return 8;
#else
    return 2;
#endif
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

#ifdef COSMO_PLANCK_15
    std::auto_ptr<PlanckLikelihood> like;
    switch(i)
    {
    case 0:
        like.reset(new PlanckLikelihood(true, false, false, false, false, false, false));
        subTestName = "LCDM_low_T";
        expected = 14.2889;
        break;
    case 1:
        like.reset(new PlanckLikelihood(true, true, false, false, false, false, false));
        subTestName = "LCDM_low_TP";
        expected = 10497.7;
        break;
    case 2:
        like.reset(new PlanckLikelihood(true, true, true, false, true, false, false));
        subTestName = "LCDM_low_TP_high_T_lite";
        expected = 10870.3;
        break;
    case 3:
        like.reset(new PlanckLikelihood(true, true, true, false, true, true, false));
        subTestName = "LCDM_low_TP_high_T_lite_lens_T";
        expected = 10879.2;
        break;
    case 4:
        like.reset(new PlanckLikelihood(true, true, true, false, true, true, true));
        subTestName = "LCDM_low_TP_high_T_lite_lens_TP";
        expected = 10886.5;
        break;
    case 5:
        like.reset(new PlanckLikelihood(true, true, true, false, false, true, true));
        subTestName = "LCDM_low_TP_high_T_lens_TP";
        expected = 16265.8;
        break;
    case 6:
        like.reset(new PlanckLikelihood(true, true, true, true, false, true, true));
        subTestName = "LCDM_low_TP_high_TP_lens_TP";
        expected = 18004.2;
        break;
    case 7:
        like.reset(new PlanckLikelihood(true, true, true, true, false, true, true, true));
        subTestName = "LCDM_Tensor_low_TP_high_TP_lens_TP";
        expected = 18011.1;
        break;
    default:
        check(false, "");
        break;
    }
    if(i < 7)
        like->setCosmoParams(paramsLCDM);
    else
        like->setCosmoParams(paramsLCDMTens);
    if(i == 5)
    {
        std::vector<double> extraParams(15);
        extraParams[0] = 100;
        extraParams[1] = -1.3;
        extraParams[2] = 0.5;
        extraParams[3] = 5;
        extraParams[4] = 200;
        extraParams[5] = 200;
        extraParams[6] = 200;
        extraParams[7] = 200;
        extraParams[8] = 5;
        extraParams[9] = 7;
        extraParams[10] = 9;
        extraParams[11] = 21;
        extraParams[12] = 80;
        extraParams[13] = 0.999;
        extraParams[14] = 0.995;
        like->setHighExtraParams(extraParams);
    }
    if(i == 6 || i == 7)
    {
        std::vector<double> extraParams(32);
        extraParams[0] = 100;
        extraParams[1] = -1.3;
        extraParams[2] = 0.5;
        extraParams[3] = 5;
        extraParams[4] = 200;
        extraParams[5] = 200;
        extraParams[6] = 200;
        extraParams[7] = 200;
        extraParams[8] = 5;
        extraParams[9] = 7;
        extraParams[10] = 9;
        extraParams[11] = 21;
        extraParams[12] = 80;
        extraParams[13] = 0.06;
        extraParams[14] = 0.05;
        extraParams[15] = 0.11;
        extraParams[16] = 0.10;
        extraParams[17] = 0.24;
        extraParams[18] = 0.72;
        extraParams[19] = -2.4;
        extraParams[20] = 0.14;
        extraParams[21] = 0.12;
        extraParams[22] = 0.30;
        extraParams[23] = 0.24;
        extraParams[24] = 0.60;
        extraParams[25] = 1.8;
        extraParams[26] = -2.4;
        extraParams[27] = 0.999;
        extraParams[28] = 0.995;
        extraParams[29] = 1;
        extraParams[30] = 1;
        extraParams[31] = 1;
        like->setHighExtraParams(extraParams);
    }
    res = like->likelihood();
#else
    PlanckLikelihood like(true, true, true, true, false, bool(i));

    if(i == 0)
        like.setCosmoParams(paramsLCDM);
    else
        like.setCosmoParams(paramsLCDMTens);
    like.setCamspecExtraParams(153, 54.9, 55.8, 4, 55.5, 4, 0.91, 0.63, 0.6, 1, 1, 0.1, 1, 0.3);

    subTestName = std::string("LCDM");
    if(i == 1)
        subTestName = std::string("LCDM+Tensor");
    res = like.likelihood();
    expected = (i == 0 ? 10038.62 : 10049.59); // calculated by running the Planck likelihood code itself on the same cl values
#endif
}

