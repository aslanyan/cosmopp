#include <vector>
#include <string>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <utils.hpp>
#include <mask_apodizer.hpp>
#include <random.hpp>
#include <histogram.hpp>
#include <likelihood.hpp>
#include <cmb.hpp>
#include <master.hpp>
#include <simulate.hpp>
#include <progress_meter.hpp>
#include <chi_squared.hpp>
#include <test_like_high.hpp>

#include <chealpix.h>
#include <healpix_base.h>
#include <alm.h>
#include <alm_healpix_tools.h>
#include <xcomplex.h>
#include <fitshandle.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>

std::string
TestLikeHigh::name() const
{
    return std::string("HIGH L LIKELIHOOD TESTER");
}

unsigned int
TestLikeHigh::numberOfSubtests() const
{
    return 1;
}

void
TestLikeHigh::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    using namespace Math;

    check(i >= 0 && i < 1, "invalid index " << i);

    const long nSide = 2048;
    const int lMaxSim = 2 * int(nSide);

    const int lMin = 31;
    const int lMax = 2000;
    const double fwhm = double(5) / 60.0;

    const int n = 1000;
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

    std::vector<double> clTT;

    output_screen1("Calculating Cl..." << std::endl);
    cmb.preInitialize(lMaxSim + 1000, false, true, false);
    cmb.initialize(paramsLCDM, true, false, true, false);
    cmb.getLensedCl(&clTT);
    output_screen1("OK" << std::endl);

    //const double nl = 0.0005;
    double nl = 0.005;

    std::vector<double> beam(lMaxSim + 1);
    Utils::readPixelWindowFunction(beam, nSide, lMaxSim, fwhm);

    Healpix_Map<double> uniformMask;
    uniformMask.SetNside(nSide, RING);
    for(unsigned long i = 0; i < uniformMask.Npix(); ++i)
        uniformMask[i] = 1;

    Healpix_Map<double>& mask = uniformMask;

    const double omegaPixel = 4 * Math::pi / mask.Npix();
    const double w = omegaPixel / nl;

    // mask out some stuff
    std::stringstream maskFileName1;
    maskFileName1 << "slow_test_files/test_highl_mask_ap_" << i << ".fits";
    Healpix_Map<double> apodizedMask;
    /*
    output_screen1("Masking out some regions..." << std::endl);
    int nRegions = 25;
    time_t seed1 = 1000000;
    Math::UniformRealGenerator thetaGen(seed1, pi / 50, pi - pi / 50), phiGen(seed1 + 1, 0, 2 * pi), angleGen(seed1 + 2, pi / 60, pi / 40);
    std::vector<double> maskTheta(nRegions), maskPhi(nRegions), maskAngle(nRegions);
    for(int i = 0; i < nRegions; ++i)
    {
        maskTheta[i] = thetaGen.generate();
        maskPhi[i] = phiGen.generate();
        maskAngle[i] = angleGen.generate();
    }

    Utils::maskRegions(mask, maskTheta, maskPhi, maskAngle);
    for(unsigned long i = 0; i < mask.Npix(); ++i)
    {
        double theta, phi;
        pix2ang_ring(nSide, i, &theta, &phi);

        if(theta < Math::pi / 2 + pi / 20 && theta > pi / 2 - pi / 20)
            mask[i] = 0;
    }
    output_screen1("OK" << std::endl);

    fitshandle outMaskHandle;
    std::stringstream maskFileName;
    maskFileName << "slow_test_files/test_highl_mask_" << i << ".fits";
    outMaskHandle.create(maskFileName.str().c_str());
    write_Healpix_map_to_fits(outMaskHandle, mask, PLANCK_FLOAT64);
    outMaskHandle.close();

    output_screen1("Apodizing the mask..." << std::endl);
    MaskApodizer ap(mask);
    const double apAngle = 30.0 / 60.0 / 180.0 * pi;
    ap.apodize(MaskApodizer::COSINE_APODIZATION, apAngle, apodizedMask);
    output_screen1("OK" << std::endl);

    fitshandle outMaskHandle1;
    outMaskHandle1.create(maskFileName1.str().c_str());
    write_Healpix_map_to_fits(outMaskHandle1, apodizedMask, PLANCK_FLOAT64);
    outMaskHandle1.close();
    */

    read_Healpix_map_from_fits(maskFileName1.str(), apodizedMask);
    //read_Healpix_map_from_fits("comb_mask_ap.fits", apodizedMask);
    //apodizedMask.swap_scheme();
    check(apodizedMask.Nside() == nSide, "");

    Healpix_Map<double> noiseWeight, noiseInvMat;
    noiseWeight.SetNside(nSide, RING);
    noiseInvMat.SetNside(nSide, RING);

    double nwAvg = 0;

    for(unsigned long i = 0; i < noiseWeight.Npix(); ++i)
    {
        double theta, phi;
        pix2ang_ring(nSide, i, &theta, &phi);

        noiseInvMat[i] = w * (1 + 0.2 * std::sin(theta) + 0.1 * std::cos(phi));
        nwAvg += 1.0 / noiseInvMat[i];

        noiseWeight[i] = noiseInvMat[i] * apodizedMask[i];
    }

    nwAvg /= noiseWeight.Npix();
    check(nwAvg > 0, "");

    nl = omegaPixel * nwAvg;

    std::vector<double> nlTT(lMaxSim + 1, nl);

    std::ofstream outCl("slow_test_files/test_like_high_cl.txt");
    for(int l = 0; l <= lMax; ++l)
    {
        const double factor = l * (l + 1) / (2 * Math::pi);
        outCl << l << '\t' << clTT[l] * factor << '\t' << nlTT[l] * factor << std::endl;
    }
    outCl.close();

    std::stringstream couplingKernelFileName;
    couplingKernelFileName << "slow_test_files/" << "test_highl_uniform_mask_cc.txt";
    std::stringstream noiseCouplingKernelFileName;
    noiseCouplingKernelFileName << "slow_test_files/" << "test_highl_noise_cc.txt";

    output_screen1("Calculating mask coupling kernel..." << std::endl);
    //Master::calculateCouplingKernel(apodizedMask, lMax + 500, couplingKernelFileName.str().c_str());
    output_screen1("OK" << std::endl);

    output_screen1("Calculating noise coupling kernel..." << std::endl);
    Master::calculateCouplingKernel(noiseWeight, lMax + 500, noiseCouplingKernelFileName.str().c_str());
    output_screen1("OK" << std::endl);

    Master master(apodizedMask, couplingKernelFileName.str().c_str(), beam, NULL, lMax + 500); 

    output_screen1("Starting simulations..." << std::endl);

    time_t seed = 10;

    StandardException exc;
    std::ofstream outLike("slow_test_files/test_highl_likes.txt");
    if(!outLike)
    {
        std::string exceptionStr = "Cannot write into file slow_test_files/test_highl_likes.txt";
        exc.set(exceptionStr);
        throw exc;
    }

    Math::Histogram<double> chi2Hist;

    Healpix_Map<double> noiseMap;
    noiseMap.SetNside(nSide, RING);

    ProgressMeter meter(n);
    for(int i = 0; i < n; ++i)
    {
        Alm<xcomplex<double> > alm, noiseAlm;
        Simulate::simulateAlm(clTT, alm, lMaxSim, seed);
        ++seed;
        //Simulate::simulateAlm(nlTT, noiseAlm, lMaxSim, seed);
        Math::GaussianGenerator noiseGen(seed, 0.0, 1.0);
        for(unsigned long j = 0; j < noiseMap.Npix(); ++j)
        {
            check(noiseInvMat[j] > 0, "");
            const double r = noiseGen.generate();
            noiseMap[j] = r * std::sqrt(1.0 / noiseInvMat[j]);
        }
        ++seed;
        noiseAlm.Set(lMaxSim, lMaxSim);
        arr<double> weight(2 * noiseMap.Nside(), 1);
        map2alm(noiseMap, noiseAlm, weight);
        for(int l = 0; l <= lMaxSim; ++l)
        {
            for(int m = 0; m <= l; ++m)
            {
                alm(l, m) += noiseAlm(l, m);
                alm(l, m) *= beam[l];
            }
        }

        Healpix_Map<double> map;
        map.SetNside(nSide, RING);
        alm2map(alm, map);

        /*
        for(unsigned long j = 0; j < map.Npix(); ++j)
            map[j] += noiseMap[j];
        */

        master.calculate(map);
        std::vector<double> clSim(lMax + 1);
        for(int l = 0; l <= lMax; ++l)
        {
            std::map<double, double> ps = master.powerSpectrum();
            check(ps.find(double(l)) != ps.end(), "");
            const double lFactor = (l == 0 ? 1 : 2.0 * Math::pi / double(l * (l + 1)));
            clSim[l] = ps[double(l)] * lFactor - nlTT[l];
        }

        LikelihoodHigh likelihood(clSim, nlTT, couplingKernelFileName.str().c_str(), lMin, lMax, noiseCouplingKernelFileName.str().c_str());
        double like;
        like = likelihood.calculate(clTT);
        outLike << like << std::endl;
        chi2Hist.addData(like);
        meter.advance();
    }
    outLike.close();
    output_screen1("OK" << std::endl);

    chi2Hist.createHistogram(chi2Hist.min(), chi2Hist.max());
    
    typedef Math::Histogram<double>::HistogramType HistogramType;
    
    const HistogramType& chi2Histogram = chi2Hist.getHistogram();
    
    output_screen1("A total of " << chi2Hist.getDataSize() << " elements in the histogram." << std::endl);
    
    output_screen1("Calculating the chi squared distribution..." << std::endl);
    std::ofstream outDist("slow_test_files/test_highl_chi2.txt");

    if(!outDist)
    {
        std::string exceptionStr = "Cannot write into file slow_test_files/test_highl_chi2.txt.";
        exc.set(exceptionStr);
        throw exc;
    }

    const double min = chi2Hist.min(), max = chi2Hist.max();
    const int dof = lMax - lMin + 1;
    Math::ChiSquared chi2(dof);
    int num = 10000;
    const double step = (max - min) / num;
    for(int i = 0; i < num; ++i)
    {
        const double x = min + step * i;
        const double y = chi2.evaluate(x);
        outDist << x << ' ' << y << std::endl;
    }
    outDist.close();
    output_screen1("OK" << std::endl);

    std::ofstream outChi2("slow_test_files/test_highl_chi2_histogram.txt");
    if(!outChi2)
    {
        std::string exceptionStr = "Cannot write into file slow_test_files/test_highl_chi2_histogram.txt.";
        exc.set(exceptionStr);
        throw exc;
    }
    
    double normalization = chi2Hist.getDataSize() * (chi2Hist.max() - chi2Hist.min()) / chi2Histogram.size();
    double myChi2 = 0;
    int myDof = 0;
    for(HistogramType::const_iterator it = chi2Histogram.begin(); it != chi2Histogram.end(); ++it)
    {
        ++myDof;
        const double x = (*it).first;
        const double y = (double)(*it).second / normalization;
        outChi2 << x << ' ' << y << std::endl;
        const double expectedY = chi2.evaluate(x);
        const double diff = y - expectedY;
    }
    
    outChi2.close();
    
    res = 1;
    expected = 1;
    subTestName = "uniform_mask";
}

