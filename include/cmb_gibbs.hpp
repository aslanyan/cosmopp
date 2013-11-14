#ifndef COSMO_CPP_CMB_GIBBS_HPP
#define COSMO_CPP_CMB_GIBBS_HPP

#include <vector>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include <healpix_map.h>
#include <alm.h>
#include <xcomplex.h>

class CMBGibbsSampler
{
public:
    typedef std::vector<std::vector<double> > GibbsSampleChain;

public:
    CMBGibbsSampler(const char* mapName, const char* noiseMapName, const char* maskName, double pixelNoise, int lMax, double fwhm, const char* startingClFileName, time_t seed = 0);
    CMBGibbsSampler(const Healpix_Map<double>& map, const Healpix_Map<double>& noiseMap, const Healpix_Map<double>& mask, double pixelNoise, int lMax, double fwhm, const std::vector<double>& startingCl, time_t seed = 0);
    ~CMBGibbsSampler();

    void step();

    const std::vector<double>& getSigmaL() const { return sigmaL_; }
    const std::vector<double>& getCl() const { return cl_; }

    static double calculateLikelihood(const std::vector<double>& cl, const std::vector<double>& sigmaL, int lMax);
    static double calculateLikelihood(const std::vector<double>& cl, const GibbsSampleChain& chain, int lMax);

    void generateChain(GibbsSampleChain& chain, int nSamples, int burnIn = 10);
    static void writeChainIntoFile(const GibbsSampleChain& chain, const char* fileName);
    static void readChainFromFile(GibbsSampleChain& chain, const char* fileName);

private:
    void calculateSigmaL();

    void generateCl();
    void generateSignal();
    void generateW();

private:
    int lMax_;
    Healpix_Map<double> map_;
    Healpix_Map<double> mask_;
    std::vector<double> cl_, sigmaL_;
    Alm<xcomplex<double> > s_;
    std::vector<double> beam_;
    double pixelNoise_;
    double a00_, a1m1_, a10_, a11_; // monopole and dipole coefficients
    Healpix_Map<double> y00_, y1m1_, y10_, y11_;

    boost::variate_generator<boost::mt19937, boost::normal_distribution<> >* generator_;
};

#endif

