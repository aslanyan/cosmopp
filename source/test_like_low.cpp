#include <vector>
#include <string>
#include <sstream>
#include <ctime>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <utils.hpp>
#include <random.hpp>
#include <histogram.hpp>
#include <c_matrix_generator.hpp>
#include <likelihood.hpp>
#include <cmb.hpp>
#include <simulate.hpp>
#include <progress_meter.hpp>
#include <chi_squared.hpp>
#include <test_like_low.hpp>

#include <chealpix.h>
#include <healpix_base.h>
#include <alm.h>
#include <alm_healpix_tools.h>
#include <xcomplex.h>
#include <fitshandle.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>

std::string
TestLikeLow::name() const
{
    return std::string("LOW L LIKELIHOOD TESTER");
}

unsigned int
TestLikeLow::numberOfSubtests() const
{
    return 1;
}

bool checkFileExists(const char* name)
{
    std::ifstream in(name);
    if(!in)
        return false;
    in.close();
    return true;
}

void
TestLikeLow::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    using namespace Math;

    check(i >= 0 && i < 1, "invalid index " << i);

    const long nSide = 16;
    const int lMaxSim = 64;

    const int lMin = 2;
    const int lMax = 30;
    const double fwhm = 10.0;

    const int n = 5000;

    const double h = 0.6704;
    const double omBH2 = 0.022032;
    const double omCH2 = 0.12038;
    const double tau = 0.0925;
    const double ns = 0.9619;
    const double as = 2.2154e-9;
    const double pivot = 0.05;

    const double pixelNoise = 1.0;

    LambdaCDMParams paramsLCDM(omBH2, omCH2, h, tau, ns, as, pivot);

    CMB cmb;

    std::vector<double> clTT;

    output_screen1("Calculating Cl..." << std::endl);
    cmb.preInitialize(lMaxSim + 1000, false, true, false);
    cmb.initialize(paramsLCDM, true, false, true, false);
    cmb.getLensedCl(&clTT);
    output_screen1("OK" << std::endl);

    std::vector<double> beam(lMaxSim + 1);
    Utils::readPixelWindowFunction(beam, nSide, lMaxSim, fwhm);

    Healpix_Map<double> uniformMask;
    uniformMask.SetNside(nSide, NEST);
    for(unsigned long i = 0; i < uniformMask.Npix(); ++i)
        uniformMask[i] = 1;

    Healpix_Map<double>& mask = uniformMask;

    // mask out some stuff
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
        pix2ang_nest(nSide, i, &theta, &phi);

        if(theta < Math::pi / 2 + pi / 20 && theta > pi / 2 - pi / 20)
            mask[i] = 0;
    }
    output_screen1("OK" << std::endl);

    fitshandle outMaskHandle;
    std::stringstream maskFileName;
    maskFileName << "slow_test_files/test_lowl_mask_" << i << ".fits";
    if(!checkFileExists(maskFileName.str().c_str()))
    {
        outMaskHandle.create(maskFileName.str().c_str());
        write_Healpix_map_to_fits(outMaskHandle, mask, PLANCK_FLOAT64);
        outMaskHandle.close();
    }

    std::vector<int> goodPixels;
    long nSideMask;
    Utils::readMask(maskFileName.str().c_str(), nSideMask, goodPixels);

    output_screen1(goodPixels.size() << " pixels left after masking." << std::endl);

    check(nSideMask == nSide, "");

    std::vector<LaVectorDouble> t(n);

    output_screen1("Starting simulations..." << std::endl);

    time_t seed = std::time(0);

    StandardException exc;

    Healpix_Map<double> noiseMap;
    noiseMap.SetNside(nSide, NEST);

    ProgressMeter meter(n);
    for(int i = 0; i < n; ++i)
    {
        Alm<xcomplex<double> > alm;
        Simulate::simulateAlm(clTT, alm, lMaxSim, seed);
        ++seed;
        Simulate::simulateWhiteNoise(noiseMap, pixelNoise, seed);
        ++seed;
        for(int l = 0; l <= lMaxSim; ++l)
        {
            for(int m = 0; m <= l; ++m)
            {
                alm(l, m) *= beam[l];
            }
        }

        Healpix_Map<double> map;
        map.SetNside(nSide, RING);
        alm2map(alm, map);
        map.swap_scheme();

        LaVectorDouble v(goodPixels.size());
        for(int j = 0; j < goodPixels.size(); ++j)
            v(j) = map[goodPixels[j]] + noiseMap[goodPixels[j]];

        t[i] = v;
        meter.advance();
    }
    output_screen1("OK" << std::endl);

    output_screen("Calculating likelihood..." << std::endl);
    std::vector<double> clCopy = clTT;
    clCopy.resize(lMax + 1);
    CMatrix* cMatrix = CMatrixGenerator::clToCMatrix(clCopy, nSide, fwhm, &goodPixels);
    CMatrix* fiducialMatrix = CMatrixGenerator::getFiducialMatrix(clTT, nSide, lMax, fwhm, &goodPixels);
    CMatrix* noiseMatrix = CMatrixGenerator::generateNoiseMatrix(nSide, pixelNoise);
    noiseMatrix->maskMatrix(goodPixels);
    LaVectorDouble foreground;
    Likelihood like(*cMatrix, *fiducialMatrix, *noiseMatrix, goodPixels, foreground);
    std::vector<std::string> mapNames(n, "test_map");
    std::vector<LikelihoodResult> results;
    like.calculateAll(t, mapNames, results);
    output_screen("OK" << std::endl);

    output_screen("Writing the output..." << std::endl);
    std::ofstream outLike("slow_test_files/test_lowl_likes.txt");
    if(!outLike)
    {
        std::string exceptionStr = "Cannot write into file slow_test_files/test_lowl_likes.txt";
        exc.set(exceptionStr);
        throw exc;
    }
    Math::Histogram<double> chi2Hist;

    for(int i = 0; i < n; ++i)
    {
        outLike << results[i].chi2 << std::endl;
        chi2Hist.addData(results[i].chi2);
    }
    outLike.close();
    output_screen("OK" << std::endl);

    chi2Hist.createHistogram(chi2Hist.min(), chi2Hist.max());
    
    typedef Math::Histogram<double>::HistogramType HistogramType;
    
    const HistogramType& chi2Histogram = chi2Hist.getHistogram();
    
    output_screen1("A total of " << chi2Hist.getDataSize() << " elements in the histogram." << std::endl);
    
    output_screen1("Calculating the chi squared distribution..." << std::endl);
    std::ofstream outDist("slow_test_files/test_lowl_chi2.txt");

    if(!outDist)
    {
        std::string exceptionStr = "Cannot write into file slow_test_files/test_lowl_chi2.txt.";
        exc.set(exceptionStr);
        throw exc;
    }

    const double min = chi2Hist.min(), max = chi2Hist.max();
    const int dof = goodPixels.size() - 1;
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

    std::ofstream outChi2("slow_test_files/test_lowl_chi2_histogram.txt");
    if(!outChi2)
    {
        std::string exceptionStr = "Cannot write into file slow_test_files/test_lowl_chi2_histogram.txt.";
        exc.set(exceptionStr);
        throw exc;
    }
    
    double normalization = chi2Hist.getDataSize() * (chi2Hist.max() - chi2Hist.min()) / chi2Histogram.size();
    const double deltaX = (chi2Hist.max() - chi2Hist.min()) / chi2Histogram.size();
    double myChi2 = 0;
    int myDof = 0;
    for(HistogramType::const_iterator it = chi2Histogram.begin(); it != chi2Histogram.end(); ++it)
    {
        const double x = (*it).first + deltaX / 2;
        const double y = (double)(*it).second / normalization;
        outChi2 << x << ' ' << y << std::endl;
        const double expectedY = chi2.evaluate(x);

        const double diff = y - expectedY;
        const double e = std::sqrt(expectedY / normalization);
        myChi2 += diff * diff / (e * e);
        ++myDof;
    }

    output_screen("Chi^2 = " << myChi2 << " per " << myDof - 1 <<  " degrees of freedom." << std::endl);
    
    outChi2.close();
    
    res = 1;
    expected = 1;

    const double chi2PerDof = myChi2 / (myDof - 1);
    const double chi2Sigma = std::sqrt(2.0 * (myDof - 1));
    if(std::abs(myChi2 - myDof +1) > 3 * chi2Sigma)
    {
        output_screen("FAIL: Not a good fit!" << std::endl);
        res = 0;
    }

    subTestName = "lowl";
}

