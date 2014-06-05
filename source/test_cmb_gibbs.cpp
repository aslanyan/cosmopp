#include <macros.hpp>
#include <cmb_gibbs.hpp>
#include <cmb.hpp>
#include <simulate.hpp>
#include <utils.hpp>
#include <test_cmb_gibbs.hpp>

#include "healpix_base.h"
#include "alm.h"
#include "alm_healpix_tools.h"
#include "xcomplex.h"
#include "healpix_map.h"
#include "rotmatrix.h"
#include "alm_powspec_tools.h"
#include "chealpix.h"

std::string
TestCMBGibbs::name() const
{
    return std::string("CMB GIBBS TESTER");
}

unsigned int
TestCMBGibbs::numberOfSubtests() const
{
    return 1;
}

void
TestCMBGibbs::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
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

    CMB cmb;

    const int lMaxSim = 2000;
    const int lMax = 64;
    const int lMaxGibbs = 30;

    std::vector<double> clTT;
    cmb.preInitialize(lMaxSim, false, true, false);
    cmb.initialize(paramsLCDM, true, false, true, false);
    cmb.getLensedCl(&clTT);

    Alm<xcomplex<double> > almSim;
    Simulate::simulateAlm(clTT, almSim, lMax);

    std::vector<double> beam(lMax + 1);
    const double fwhm = 10;
    const long nSide = 32;
    const double pixelNoise = 1.0;

    Utils::readPixelWindowFunction(beam, nSide, lMax, fwhm);

    for(int l = 0; l <= lMax; ++l)
    {
        for(int m = 0; m <= l; ++m)
        {
            almSim(l, m) *= beam[l];
        }
    }

    Healpix_Map<double> map, noiseMap, uniformMask;
    map.SetNside(nSide, RING);
    alm2map(almSim, map);
    map.swap_scheme();

    noiseMap.SetNside(nSide, RING);
    uniformMask.SetNside(nSide, RING);
    for(long i = 0; i < noiseMap.Npix(); ++i)
    {
        noiseMap[i] = pixelNoise;
        uniformMask[i] = 1;
    }
    noiseMap.swap_scheme();
    uniformMask.swap_scheme();

    std::vector<double> clStarting(lMaxGibbs + 1, 1000.0);
    for(int l = 1; l <= lMaxGibbs; ++l)
        clStarting[l] *= (2 * Math::pi / (l * (l + 1)));
    clStarting[0] = clStarting[1];

    CMBGibbsSampler gibbs(map, noiseMap, uniformMask, pixelNoise, lMaxGibbs, fwhm, clStarting);
    CMBGibbsSampler::GibbsSampleChain chain;

    gibbs.generateChain(chain, 100);

    switch(i)
    {
    case 0:
        subTestName = "uniform_mask";
        res = 1;
        expected = 1;

        break;
    default:
        check(false, "");
        break;
    }
}

