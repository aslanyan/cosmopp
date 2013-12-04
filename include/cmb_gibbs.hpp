#ifndef COSMO_CPP_CMB_GIBBS_HPP
#define COSMO_CPP_CMB_GIBBS_HPP

#include <vector>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include <healpix_map.h>
#include <alm.h>
#include <xcomplex.h>

/// A CMB Gibbs Sampler.
class CMBGibbsSampler
{
public:
    /// The Gibbs chain type.
    typedef std::vector<std::vector<double> > GibbsSampleChain;

public:
    /// Constructor.
    /// \param mapName The name of the file containing the map to be sampled.
    /// \param noiseMapName The name of the file containg the noise map to be added to the data map (a 0 map should be passed here if the noise is already in the data).
    /// \param maskName The name of the mask file.
    /// \param pixelNoise The uniform noise per pixel.
    /// \param lMax The maximum value of l for sampling.
    /// \param fwhm The full width at half maximum of the beam function of the map (in degrees).
    /// \param startingClFileName The name of the file containing the C_l values for starting the sampling. Each line of the file should contain one C_l value, starting from l = 0.
    /// \param seed A random seed for scanning. If the value is 0 (which it is by default) the seed will be chosen from the current time.
    CMBGibbsSampler(const char* mapName, const char* noiseMapName, const char* maskName, double pixelNoise, int lMax, double fwhm, const char* startingClFileName, time_t seed = 0);

    /// Constructor.
    /// \param map The map to be sampled.
    /// \param noiseMap The map to be added to the data map (a 0 map should be passed here if the noise is already in the data).
    /// \param mask The mask.
    /// \param pixelNoise The uniform noise per pixel.
    /// \param lMax The maximum value of l for sampling.
    /// \param fwhm The full width at half maximum of the beam function of the map (in degrees).
    /// \param startingCl A vector containing the C_l values for starting the sampling, starting from l = 0.
    /// \param seed A random seed for scanning. If the value is 0 (which it is by default) the seed will be chosen from the current time.
    CMBGibbsSampler(const Healpix_Map<double>& map, const Healpix_Map<double>& noiseMap, const Healpix_Map<double>& mask, double pixelNoise, int lMax, double fwhm, const std::vector<double>& startingCl, time_t seed = 0);

    /// Destructor.
    ~CMBGibbsSampler();

    /// Take one step in the sampling.
    void step();

    /// Get the current sample for sigma_l = sum_m |a_lm|^2.
    /// \return A vector containg the sigma_l values, starting from l = 0.
    const std::vector<double>& getSigmaL() const { return sigmaL_; }

    /// Get the current C_l sample.
    /// \return A vector containing the C_l values, starting from l = 0.
    const std::vector<double>& getCl() const { return cl_; }

    /// Calculate the likelihood of given C_l values given a sigma_l sample.
    /// \param cl A vector containg the C_l values for which the likelihood must be calculated, starting from l = 0.
    /// \param sigmaL A vector containing the sigma_l sample, starting from l = 0.
    /// \param lMax The maximum value of l to be used in the calculation. Both cl and sigmaL must have sizes >= lMax + 1.
    /// \return -2ln(likelihood).
    static double calculateLikelihood(const std::vector<double>& cl, const std::vector<double>& sigmaL, int lMax);

    /// Calculate the likelihood of given C_l values given a Gibbs chain.
    /// \param cl A vector containg the C_l values for which the likelihood must be calculated, starting from l = 0.
    /// \param chain The Gibbs chain to be used for the likelihood calculation, can be generated using generateChain.
    /// \param lMax The maximum value of l to be used in the calculation. The chain elements and cl must have sizes >= lMax + 1.
    /// \return -2ln(likelihood).
    static double calculateLikelihood(const std::vector<double>& cl, const GibbsSampleChain& chain, int lMax);

    /// Generate a Gibbs chain. The chain will be generated from the current state of the sampler, i.e. samples will be generated and written in the chain starting from the current state.
    /// \param chain The Gibbs chain will be written here.
    /// \param nSamples The number of samples that the chain must contain.
    /// \param burnIn The number of samples to ignore before starting to write the chain.
    void generateChain(GibbsSampleChain& chain, int nSamples, int burnIn = 10);

    /// Write a given Gibbs chain into a file.
    /// \param chain The chain to be written in a file.
    /// \param fileName The name of the file to contain the chain.
    static void writeChainIntoFile(const GibbsSampleChain& chain, const char* fileName);

    /// Read a Gibbs chain from a file.
    /// \param chain The chain will be written here.
    /// \param fileName The name of the file to read from.
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

