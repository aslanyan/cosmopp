#ifndef COSMO_PP_LIKELIHOOD_HPP
#define COSMO_PP_LIKELIHOOD_HPP

#include <vector>
#include <string>

#include <complex_types.hpp>
#include <c_matrix.hpp>
#include <whole_matrix.hpp>
#include <matrix.hpp>

#include "alm.h"
#include "xcomplex.h"

/// A result of the Likelihood calculation.

/// This is a simple struct for the likelihood calculation result.
struct LikelihoodResult
{
    std::string mapName; ///< The name of the map.
    double logDet; ///< Logarithm of the determinant of the covariance matrix.
    double chi2; ///< Chi squared.
    double like; ///< -2Log(likelihood), this is simply the sum of logDet and chi2.
    
    /// Compares two results.
    inline bool operator < (const LikelihoodResult& other) const { return like < other.like; }
    
    /// Copy another result.
    LikelihoodResult& operator = (const LikelihoodResult& other) { mapName = other.mapName; logDet = other.logDet; chi2 = other.chi2; like = other.like; return *this; }
};

/// Likelihood calculation class for temperature maps.

/// This class is used to calculate the likelihood for temperature maps. Foreground marginalization is included.
/// The calculation is in pixel space, can be used only for low-l.
class Likelihood
{
public:
    /// Constructor.
    
    /// Constructs a likelihood calculator.
    /// \param cMatrix The covariance matrix.
    /// \param fiducialMatrix The fiducial covariance matrix.
    /// \param noiseMatrix The noise covariance matrix.
    /// \param maskFileName The file name of the mask to be used.
    /// \param foregroundFileName The file name of a foreground map (fits format) to marginalize over. NULL value implies no foreground marginalization.
    Likelihood(const CMatrix& cMatrix, const CMatrix& fiducialMatrix, const CMatrix& noiseMatrix, const char* maskFileName, const char* foregroundFileName = NULL);
    
    
    /// Constructor.
    
    /// Constructs a likelihood calculator.
    /// \param cMatrix The covariance matrix.
    /// \param fiducialMatrix The fiducial covariance matrix.
    /// \param noiseMatrix The noise covariance matrix.
    /// \param goodPixels A vector containing the indices of the unmasked pixels.
    /// \param foreground A vector containing the foreground map. This can be read by the function readForeground. If it's empty foreground marginalization is not done.
    Likelihood(const CMatrix& cMatrix, const CMatrix& fiducialMatrix, const CMatrix& noiseMatrix, const std::vector<int>& goodPixels, const std::vector<double>& foreground);
    
    /// Calculate likelihood for a single map.
    
    /// This function calculates the likelihood function for a single map.
    /// \param mapName The name of the map file (fits).
    /// \param noiseMapName The name of the file that contains noise to be added to the map (useful for simulations etc.).
    /// \param chi2 The chi squared value to be returned.
    /// \param logDet The logarithm of the determinant of the covariance matrix to be returned.
    /// \return -2Log(likelihood) which is the sum of chi2 and logDet.
    double calculate(const char* mapName, const char* noiseMapName, double& chi2, double& logDet) const;
    
    /// Calculate likelihood for many maps.
    
    /// This function calculates the likelihood for many maps.
    /// \param inputListName The name of a file containing all the map names. The format is, the first line contains the number of maps, each following line contains the name of the map followed by the name of the noise map.
    /// \param results A vector to return the calculated results in. The new results are added to what exists in the vector.
    void calculateAll(const char* inputListName, std::vector<LikelihoodResult>& results) const;
    
    /// Calculate likelihood for many maps.
    
    /// This function calculates the likelihood for many maps.
    /// \param t A vector of temperature maps (noise added). This can be read by the function readInput.
    /// \param mapNames A vector with the names of the maps. This can be read by the function readInput.
    /// \param results A vector to return the calculated results in. The new results are added to what exists in the vector.
    void calculateAll(const std::vector<std::vector<double> >& t, const std::vector<std::string>& mapNames, std::vector<LikelihoodResult>& results) const;
    
    /// Calculate likelihood for a single map.
    
    /// This function calculates the likelihood function for a single map.
    /// \param t A vector containing the map (noise added). This can be read by the function readMapAndNoise.
    /// \param chi2 The chi squared value to be returned.
    /// \param logDet The logarithm of the determinant of the covariance matrix to be returned.
    /// \return -2Log(likelihood) which is the sum of chi2 and logDet.
    double calculate(const std::vector<double>& t, double& chi2, double& logDet) const;
    
    /// Read input from many maps.
    
    /// This function reads many maps and noise maps and adds then together.
    /// \param inputListName The name of a file containing all the map names. The format is, the first line contains the number of maps, each following line contains the name of the map followed by the name of the noise map.
    /// \param goodPixels A vector containing the indices of the unmasked pixels.
    /// \param t A vector where the maps are returned in.
    /// \param mapNames A vector where the map names are returned in.
    static void readInput(const char* inputListName, const std::vector<int>& goodPixels, std::vector<std::vector<double> >& t, std::vector<std::string>& mapNames);
    
    /// Read a single map and noise map.
    
    /// This function reads a map and a noise map and adds them together.
    /// \param mapName The name of the map file (fits format).
    /// \param noiseMapName The name of the noise map file (fits format).
    /// \param goodPixels A vector containing the indices of unmasked pixels.
    /// \param nSide NSide of the map to be returned.
    /// \param t The map (with added noise) to be returned.
    static void readMapAndNoise(const char* mapName, const char* noiseMapName, const std::vector<int>& goodPixels, long& nSide, std::vector<double>& t);
    
    /// Read the foreground map.
    
    /// This function reads the foreground map.
    /// \param foregroundFileName The map file name (fits format).
    /// \param goodPixels A vector containing the unmasked pixels.
    /// \param nSide NSide of the foreground map to be returned.
    /// \param f The foreground map to be returned.
    static void readForeground(const char* foregroundFileName, const std::vector<int>& goodPixels, long& nSide, std::vector<double>& f);
    
private:
    void construct(const CMatrix& cMatrix, const CMatrix& fiducialMatrix, const CMatrix& noiseMatrix, const char* maskFileName, const char* foregroundFileName);
    void construct(const CMatrix& cMatrix, const CMatrix& fiducialMatrix, const CMatrix& noiseMatrix, const std::vector<int>& goodPixels, const std::vector<double>& foreground);
    double vmv(int n, const std::vector<double>& a, const Math::SymmetricMatrix<double>& matrix, const std::vector<double>& b) const; // transpose(a) * matrix * b
    
private:
    Math::SymmetricMatrix<double> cInv_;
    double logDet_;
    std::vector<double> f_;
    std::vector<int> goodPixels_;
};

/// Likelihood calculation class for polarization maps.

/// This class is used to calculate likelihood for polarization maps.
class LikelihoodPolarization
{
public:
    typedef Alm<xcomplex<double> > AlmType;
    
public:
    /// Constructor.
    
    /// Constructs a polarization likelihood calculator.
    /// The consructor needs to read the inverse noise matrix from the text file n_inv.txt. The file must have double the number of the pixels in rows and columns, the first half for Q, second half for U.
    /// \param cMatrix The covariance matrix which must include both Q and U maps. This matrix is already rotated.
    /// \param nSide Healpix NSide for Q and U maps.
    /// \param goodPixels A vector containing the unmasked pixels.
    /// \param lMax Maximum l to include in the calculation.
    /// \param etttInverse The whole matrix for ET(TT)^(-1), this must be UNROTATED.
    /// \param phi Euler angle phi of the rotation. This is needed because etttInverse is not rotated.
    /// \param theta Euler angle theta of the rotation. This is needed because etttInverse in not rotated.
    /// \param psi Euler angle psi of the rotation. This is needed because etttInverse is not rotated.
    LikelihoodPolarization(const CMatrix& cMatrix, long nSide, const std::vector<int>& goodPixels, int lMax, const WholeMatrix& etttInverse, double phi = 0, double theta = 0, double psi = 0);
    
    /// Calculates the likelihood.
    
    /// This function is used to calculate the likelihood from maps that have already been read.
    /// \param v A vector that contains the Q and U maps. Can be read using the function readMaps. IMPORTANT: The Q and U maps are not the original maps, they need to have the inverse noise matrix applied to them, i.e. they are N^(-1)*(Q, U).
    /// \param alm The Alm-s for TT, can be read using the function readMaps.
    /// \param chi2 The chi squared value to be returned.
    /// \param logDet The logarithm of the determinant of the covariance matrix to be returned.
    /// \return -2Log(likelihood) which is the sum of chi2 and logDet.
    double calculate(const std::vector<double>& v, const AlmType& alm, double& chi2, double& logDet) const;
    
    /// Calculates the likelihood.
    
    /// This function is used to calculate the likelihood from Q and U maps.
    /// \param qMap The Q map filename (fits format). IMPORTANT: The Q and U maps are not the original maps, they need to have the inverse noise matrix applied to them, i.e. they are N^(-1)*(Q, U).
    /// \param uMap The U map filename (fits format). IMPORTANT: The Q and U maps are not the original maps, they need to have the inverse noise matrix applied to them, i.e. they are N^(-1)*(Q, U).
    /// \param almTTFileName The filename where the TT Alm-s are written (text format). Complex numbers written in format (re, im) starting from l=0, m=0. It should include only m>=0 values.
    /// \param chi2 The chi squared value to be returned.
    /// \param logDet The logarithm of the determinant of the covariance matrix to be returned.
    /// \return -2Log(likelihood) which is the sum of chi2 and logDet.
    double calculate(const char* qMap, const char* uMap, const char* almTTFileName, double& chi2, double& logDet) const;
    
    /// Calculates the likelihood for a set of maps.
    
    /// This function is used to calculate the likelihood for many maps at once.
    /// \param v A vector of vectors containing Q and U maps. Each vector can be read using the function readMaps.
    /// \param alm A vector of TT Alm-s. Each Alm can be read using the function readMaps.
    /// \param mapNames A vector containig the names of maps. These names are used to write the results in the output.
    /// \param results A vector of results to be returned.
    void calculateAll(const std::vector<std::vector<double> >& v, const std::vector<AlmType>& alm, const std::vector<std::string>& mapNames, std::vector<LikelihoodResult>& results) const;
    
    /// Get an inverse noise matrix element.
    
    /// This function is used to get an element of the inverse noise matrix (Q and U combined, mask applied).
    /// \param i First index.
    /// \param j Second index.
    /// \return The inverse noise matrix element.
    double getNInv(int i, int j) const;
    
    /// Combines whole matrics.
    
    /// This function is used to produced combined (ET(TT)^(-1)TE) and ET(TT)^(-1) matrices from TT, TE, and EE matrices.
    /// \param tt Input TT whole matrix.
    /// \param te Input TE whole matrix.
    /// \param ee Input EE whole matrix.
    /// \param combined Output combined whole matrix.
    /// \param etttInverse Output ET(TT)^(-1) whole matrix.
    static void combineWholeMatrices(const WholeMatrix& tt, const WholeMatrix& te, const WholeMatrix& ee, WholeMatrix& combined, WholeMatrix& etttInverse);
    
    /// Read polarization maps and TT Alm-s.
    
    /// This function reads polarization maps and TT Alm-s.
    /// \param qMapName The Q map filename (fits format).
    /// \param uMapName The U map filename (fits format).
    /// \param almTTFileName The filename where the TT Alm-s are written (text format). Complex numbers written in format (re, im) starting from l=0, m=0. It should include only m>=0 values.
    /// \param lMax Maximum value of l to be used in likelihood calculation.
    /// \param goodPixels A vector containing the unmasked pixels.
    /// \param v A vector that will contain the Q and U maps upon return.
    /// \param alm TT Alm-s to be returned.
    static void readMaps(const char* qMapName, const char* uMapName, const char* almTTFileName, int lMax, const std::vector<int>& goodPixels, std::vector<double>& v, AlmType& alm);
    
    /// Read input from many maps.
    
    /// This function is a generalization of Likelihood::readInput. Reads many temperature and polarization maps.
    /// \param inputListName The name of a file containing all the map names. The format is, the first line contains the number of maps, each following line contains the name of the temperature map followed by the name of the noise map, TT alm file, Q map, and U map.
    /// \param goodPixels A vector containing the indices of unmasked pixels for temperature.
    /// \param t A vector that will contain temperature maps upon return.
    /// \param mapNames A vector that will contain temperature map names upon return.
    /// \param goodPixelsPol A vector containing the indices of unmasked pixels for Q and U maps.
    /// \param v A vector that will contain Q and U maps upon return.
    /// \param alm A vector that will contain TT Alm-s upon return.
    /// \param lMaxPol Maximum value of l for polarization analysis.
    static void readInput(const char* inputListName, const std::vector<int>& goodPixels, std::vector<std::vector<double> >& t, std::vector<std::string>& mapNames, const std::vector<int>& goodPixelsPol, std::vector<std::vector<double> >& v, std::vector<AlmType>& alm, int lMaxPol);
    
private:
    Math::SymmetricMatrix<double> cInv_;
    double logDet_;
    std::vector<int> goodPixels_;
    long nSide_;
    int lMax_;
    
    double phi_;
    double theta_;
    double psi_;
    
    Math::SymmetricMatrix<double> nInv_;
    
    WholeMatrix etttInverse_;
    std::vector<double> beam_;
};

/// Likelihood calculation for temperature maps.

/// This class calculates the likelihood function for temperature maps in harmonic space. Is only accurate for high l values.
class LikelihoodHigh
{
public:
    /// Constructor.
    /// \param cl Data Cl-s in muK^2. The index of the vector is l.
    /// \param nl Noise Cl-s in muK^2. The index of the vector is l.
    /// \param couplingKernelFileName The name of the file containig the coupling kernel for the mask. This is produced when running Master.
    /// \param lMin The minimum value of l to include in the calculation.
    /// \param lMax The maximum value of l to include in the calculation.
    LikelihoodHigh(const std::vector<double>& cl, const std::vector<double>& nl, const char* couplingKernelFileName, int lMin, int lMax, const char* noiseCouplingKernelFileName = 0);

    /// Constructor.
    /// \param dataClFileName The name of the file containing data Cl-s. The units should be muK^2. The format is each row should contain l followed by Cl.
    /// \param noiseClFileName The name of the file containing noise Cl-s. The units should be muK^2. The format is each row should contain l followed by Cl.
    /// \param couplingKernelFileName The name of the file containig the coupling kernel for the mask. This is produced when running Master.
    /// \param lMin The minimum value of l to include in the calculation.
    /// \param lMax The maximum value of l to include in the calculation.
    LikelihoodHigh(const char* dataClFileName, const char* noiseClFileName, const char* couplingKernelFileName, int lMin, int lMax, const char* noiseCouplingKernelFileName = 0);

    /// Calculate the likelihood for given Cl-s.
    /// \param cl Model Cl-s in muK^2. The index of the vector is l.
    /// \return -2ln(likelihood).
    double calculate(const std::vector<double>& cl) const;

    /// Calculate the likelihood for given Cl-s.
    /// \param clFileName The name of the file containing model Cl-s in muK^2. Each row should contain just the Cl, starting from l = 0.
    /// \return -2ln(likelihood).
    double calculate(const char* clFileName) const;

private:
    void readCouplingKernel(const char* couplingKernelFileName, std::vector<std::vector<double> > & coupling);

private:
    std::vector<double> cl_, nl_;
    std::vector<std::vector<double> > coupling_;
    std::vector<std::vector<double> > noiseCoupling_;
    const int lMin_, lMax_;
    const double offset_;
};

#endif
