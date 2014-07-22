#include <fstream>
#include <cmath>
#include <ctime>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <cmb_gibbs.hpp>
#include <cmb.hpp>
#include <simulate.hpp>
#include <utils.hpp>
#include <random.hpp>
#include <markov_chain.hpp>
#include <numerics.hpp>
#include <progress_meter.hpp>
#include <test_cmb_gibbs.hpp>

#include <chealpix.h>
#include <healpix_base.h>
#include <alm.h>
#include <alm_healpix_tools.h>
#include <xcomplex.h>
#include <fitshandle.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>

std::string
TestCMBGibbs::name() const
{
    return std::string("CMB GIBBS TESTER");
}

unsigned int
TestCMBGibbs::numberOfSubtests() const
{
    return 5;
}

void
TestCMBGibbs::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
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

    CMB cmb;

    const int lMaxSim = 2000;
    const int lMax = (i == 0 ? 64 : 128);
    const int lMaxGibbs = (i == 0 ? 30 : 50);

    std::vector<double> clTT;
    cmb.preInitialize(lMaxSim, false, true, false);
    cmb.initialize(paramsLCDM, true, false, true, false);
    cmb.getLensedCl(&clTT);


    const int seed = std::time(0);
    Math::UniformRealGenerator gen(seed, 0, 1000000);

    Alm<xcomplex<double> > almSim;
    Simulate::simulateAlm(clTT, almSim, lMax, int(gen.generate()));

    std::vector<double> beam(lMax + 1);
    double fwhm = (i == 0 ? 10 : 5);
    long nSide = (i == 0 ? 16 : 32);
    double pixelNoise = (i == 0 ? 1.0 : 2.0);
    if(i == 2)
        pixelNoise = 10.0;    

    if(i == 3)
        pixelNoise = 50.0;

    double comparisonScale = 0.25;
    int comparisonLMax = lMaxGibbs;
    if(i == 2)
        comparisonScale = 0.5;
    if(i == 3)
    {
        comparisonScale = 10.0;
        comparisonLMax = 42;
    }
    if(i >= 4)
    {
        comparisonScale = 5.0;
    }

    Utils::readPixelWindowFunction(beam, nSide, lMax, fwhm);

    std::vector<double> clOriginal(lMax + 1, 0);

    StandardException exc;
    std::stringstream originalClFileName;
    originalClFileName << "slow_test_files/cmb_gibbs_" << i << "_cl_original.txt";
    std::ofstream out(originalClFileName.str().c_str());
    if(!out)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into file " << originalClFileName.str() << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    for(int l = 0; l <= lMax; ++l)
    {
        for(int m = 0; m <= l; ++m)
        {
            clOriginal[l] += (m ? 2 : 1) * almSim(l, m).norm();
            almSim(l, m) *= beam[l];
        }

        clOriginal[l] /= (2 * l + 1);

        out << l << '\t' << clOriginal[l] * l * (l + 1) / (2 * Math::pi) << '\t' << clTT[l] * l * (l + 1) / (2 * Math::pi) << std::endl;
    }
    out.close();

    Healpix_Map<double> map, noiseMap, uniformMask, mask;
    map.SetNside(nSide, RING);
    alm2map(almSim, map);

    if(i == 4)
        read_Healpix_map_from_fits(std::string("slow_test_files/mask1.fits"), mask);

    noiseMap.SetNside(nSide, RING);
    uniformMask.SetNside(nSide, RING);
    for(long i = 0; i < noiseMap.Npix(); ++i)
    {
        uniformMask[i] = 1;
    }
    Simulate::simulateWhiteNoise(noiseMap, pixelNoise, int(gen.generate()));

    std::vector<double> clStarting(lMaxGibbs + 1, 1000.0);
    for(int l = 1; l <= lMaxGibbs; ++l)
        clStarting[l] *= (2 * Math::pi / (l * (l + 1)));
    clStarting[0] = clStarting[1];

    CMBGibbsSampler gibbs(map, noiseMap, (i < 4 ? uniformMask : mask), pixelNoise, lMaxGibbs, fwhm, clStarting, int(gen.generate()));

    std::stringstream chainFileName;
    chainFileName << "slow_test_files/cmb_gibbs_" << i << "_chain.txt";
    out.open(chainFileName.str().c_str());
    if(!out)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into file " << chainFileName.str() << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    const int burnin = 10;
    const int samples = 1000;

    std::vector<Posterior1D> clPosteriors;
    output_screen("Generating the Gibbs chain..." << std::endl);
    ProgressMeter meter(samples + burnin);
    for(int i = 0; i < samples + burnin; ++i)
    {
        gibbs.step();
        out << 1 << '\t' << 1;
        for(int l = 0; l <= lMaxGibbs; ++l)
            out << '\t' << gibbs.getCl()[l] * l * (l + 1) / (2 * Math::pi);
        out << std::endl;
        meter.advance();
    }
    out.close();

    MarkovChain chain(chainFileName.str().c_str(), burnin);

    std::stringstream statsFileName;
    statsFileName << "slow_test_files/cmb_gibbs_" << i << "_stats.txt";
    out.open(statsFileName.str().c_str());
    if(!out)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into file " << statsFileName.str() << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    subTestName = "uniform_mask_low_l";
    if(i == 1)
        subTestName = "uniform_mask_high_l";
    if(i == 2)
        subTestName = "uniform_mask_high_l_high_noise";
    if(i == 3)
        subTestName = "uniform_mask_high_l_very_high_noise";
    if(i == 4)
        subTestName = "custom_mask_high_l";
    res = 1;
    expected = 1;

    for(int l = 0; l <= lMaxGibbs; ++l)
    {
        double smoothScale = 0;
        if(l >= 2)
            smoothScale = 200;
        if(l >= 4)
            smoothScale = 100;
        if(l >= 10)
            smoothScale = 0;

        output_screen1("Generating the posterior distribution for l = " << l << "..." << std::endl);
        Posterior1D* post = chain.posterior(l, (l > 1 ? Posterior1D::GAUSSIAN_SMOOTHING : Posterior1D::SPLINE_SMOOTHING), smoothScale);
        output_screen1("OK" << std::endl);
        double l1, l2, u1, u2;
        post->get1SigmaTwoSided(l1, u1);
        post->get2SigmaTwoSided(l2, u2);
        const double peak = post->peak();
        out << l << '\t' << post->mean() << '\t' << peak << '\t' << post->median() << '\t' << post->median() - l1 << '\t' << post->median() + u1 << '\t' << post->median() - l2 << '\t' << post->median() + u2 << std::endl;

        const double clOrig = clOriginal[l] * l * (l + 1) / (2 * Math::pi);

        if(l <= comparisonLMax && !Math::areEqual(clOrig, peak, comparisonScale))
        {
            output_screen("FAIL: for l = " << l << " original Cl = " << clOrig << ", resulting Cl from the sample = " << peak << ". Too far away!" << std::endl);
            res = 0;
        }

        if(l >= 2)
        {
            std::stringstream postFileName;
            postFileName << "slow_test_files/cmb_gibbs_" << i << "_c_" << l << ".txt";
            std::ofstream outPost(postFileName.str().c_str());
            if(!outPost)
            {
                std::stringstream exceptionStr;
                exceptionStr << "Cannot write into file " << postFileName.str() << ".";
                exc.set(exceptionStr.str());
                throw exc;
            }
            const int nPoints = 100;
            double min = post->min(), max = post->max();
            if(max > 5000)
                max = 5000;
            if(min >= 5000)
            {
                outPost.close();
                delete post;
                continue;
            }
            const double delta = (max - min) / nPoints;
            for(int i = 0; i <= nPoints; ++i)
            {
                double t = min + delta * i;
                if(i == nPoints)
                    t = max;
                outPost << t << '\t' << post->evaluate(t) << std::endl;
            }
            outPost.close();
        }
        delete post;    
    }
    out.close();
}

