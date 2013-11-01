#ifndef COSMO_CPP_C_MATRIX_GENERATOR_HPP
#define COSMO_CPP_C_MATRIX_GENERATOR_HPP

#include <vector>

#include "c_matrix.hpp"
#include "whole_matrix.hpp"

/// A container for Legendre Polynomials calculated between pixels. Can be used in the CMatrixGenerator class to speed up calculations.
class LegendrePolynomialContainer
{
public:
    /// Constructor.
    /// \param lMax Maximum value of l to store Legendre polynomials for.
    /// \param nSide NSide of the pixelization.
    /// \param goodPixels A pointer to a vector containig the indices of unmasked pixels, NULL to use them all. Legendre polynomials are calculated between all of the unmasked pixels.
    LegendrePolynomialContainer(int lMax, long nSide, const std::vector<int>* goodPixels = NULL);
    
    /// Constructor that reads from a file.
    /// \param fileName The name of the file containing the Legendre polynomials.
    LegendrePolynomialContainer(const char* fileName);
    
    /// Retreives the value of Legendre polynomial l applied to cos(theta) where theta is the angle between pixel i and j.
    /// \param l The l value of the polynomial.
    /// \param i The index of the first pixel (indices are the same as in goodPixels given in the constructor).
    /// \param j The index of the second pixel (indices are the same as in goodPixels given in the constructor).
    double value(int l, int j, int i) const;
    
    /// Saves the contents into a file.
    /// \param fileName The name of the file to contain the Legendre polynomials.
    void writeIntoFile(const char* fileName) const;
    
private:
    //first index is l, second index is j, third index is i, each element is P_l(n_j dot n_i) * B_l * B_l, B_l is the beam function
    std::vector<std::vector<std::vector<double> > > data_;
};

/// Converts covariance matrices from l-m space to pixel space.

/// This class provides functions for converting covariance matrices from l-m space to pixel space.
class CMatrixGenerator
{
public:
    /// Creates a covariance matrix from C_l values read from a text file.
    
    /// This function creates a CMatrix from given values of C_l.
    /// The input C_l is in units of muK, the CMatrix generated has units mK.
    /// \param cl A vector with values of C_l, the index is l. The format is the same as CAMB output (see camb.info/readme.html). l_max = size of cl - 1.
    /// \param nSide NSide of the output matrix.
    /// \param fwhm The full width at half maximum of the gaussian beam.
    /// \param goodPixels A pointer to a vector containing the indices of unmasked pixels, NULL to use them all.
    /// \param lp A pointer to a container of Legendre polynomials, NULL makes it calculate them by itself.
    /// \return A pointer to the generated covariance matrix. The units are mK. It must be deleted after using.
    static CMatrix* clToCMatrix(const std::vector<double>& cl, long nSide, double fwhm, const std::vector<int>* goodPixels = NULL, const LegendrePolynomialContainer* lp = NULL);
    
    /// Creates a covariance matrix from C_l values read from a text file.
    
    /// This function creates a CMatrix from values of C_l read from a file.
    /// The input file has C_l in units of muK, the CMatrix generated has units mK.
    /// \param clFileName The name of the file to read the C_l values from. The format is the same as CAMB output (see camb.info/readme.html).
    /// \param nSide NSide of the output matrix.
    /// \param lMax Values of l up to this (starting from 2) are included in C_l.
    /// \param fwhm The full width at half maximum of the gaussian beam.
    /// \param goodPixels A pointer to a vector containing the indices of unmasked pixels, NULL to use them all.
    /// \param lp A pointer to a container of Legendre polynomials, NULL makes it calculate them by itself.
    /// \return A pointer to the generated covariance matrix. The units are mK. It must be deleted after using.
    static CMatrix* clToCMatrix(const char* clFileName, long nSide, int lMax, double fwhm, const std::vector<int>* goodPixels = NULL, const LegendrePolynomialContainer* lp = NULL);
    
    /// Generates TT, TE, and EE whole matrices from C_l values read from a text file.
    
    /// This function reads a C_l file and writes the corresponding values into TT, TE, and EE whole matrices.
    /// The l range is determined from the whole matrices upon input (they must have the same lMin and lMax).
    /// \param clFileName The name of the file to read the C_l values from. The format is the same as CAMB output (see camb.info/readme.html).
    /// \param tt The TT whole matrix where the corresponding c_l values will be written upon return. Off diagonal elements will be 0.
    /// \param te The TE whole matrix where the corresponding c_l values will be written upon return. Off diagonal elements will be 0.
    /// \param ee The EE whole matrix where the corresponding c_l values will be written upon return. Off diagonal elements will be 0.
    static void clToWholeMatrix(const char* clFileName, WholeMatrix& tt, WholeMatrix& te, WholeMatrix& ee);
    
    /// Convert temperature-temperature WholeMatrix to CMatrix.
    
    /// This function converts a given WholeMatrix to CMatrix, also including a beam function. The units are converted from muK to mK.
    /// It also rotates passively the coordinate frame by given three Euler angles. 
    /// The rotation is done as follows. First rotate counterclockwise by angle phi, then rotate counterclockwise around the new x axis by angle theta, 
    /// then rotate counterclockwise around the new z axis by angle psi.
    /// \param wholeMatrix The WholeMatrix to be converted. The units are muK.
    /// \param nSide The NSide of the pixel space to convert to.
    /// \param fwhm Full width at half maximum of the gaussian beam, in degrees.
    /// \param phi Euler angle phi.
    /// \param theta Euler angle theta.
    /// \param psi Euler angle psi.
    /// \return A pointer to the generated CMatrix. The units are mK. It must be deleted after using.
    static CMatrix* wholeMatrixToCMatrix(const WholeMatrix& wholeMatrix, long nSide, double fwhm, double phi = 0, double theta = 0, double psi = 0);
    
    /// Convert polarization E-E WholeMatrix to CMatrix.
    
    /// This function converts a given WholeMatrix to CMatrix, also including a beam function. The units are converted from muK to mK. 
    /// It also rotates passively the coordinate frame by given three Euler angles. 
    /// The rotation is done as follows. First rotate counterclockwise by angle phi, then rotate counterclockwise around the new x axis by angle theta, 
    /// then rotate counterclockwise around the new z axis by angle psi.
    /// \param ee The WholeMatrix E-E to be converted. The units are muK.
    /// \param nSide The NSide of the pixel space to convert to.
    /// \param fwhm Full width at half maximum of the gaussian beam, in degrees.
    /// \param phi Euler angle phi.
    /// \param theta Euler angle theta.
    /// \param psi Euler angle psi.
    /// \return A pointer to the generated CMatrix. The units are mK. This CMatrix has twice the number of pixels, the first half correspond to Q, the second half to U. It must be deleted after using.
    static CMatrix* polarizationEEWholeMatrixToCMatrix(const WholeMatrix& ee, long nSide, double fwhm, double phi = 0, double theta = 0, double psi = 0);
    
    /// Creates a fiducial matrix from C_l values read from a text file.
    
    /// A fiducial matrix serves two purposes. Firstly, it marginalizes over the monopole and the dipole term by having large variance terms for those. Secondly, it includes terms above a given lMax
    /// up to 4*NSide. If one is testing a model with off-diagonal elements up to a given lMax, the fiducial matrix needs to be added which will include the standard terms above that lMax.
    /// The input file has C_l in units of muk, the CMatrix generated has units mK.
    /// \param clFileName The name of the file to read the C_l values from. The format is the same as CAMB output (see camb.info/readme.html).
    /// \param nSide NSide of the output matrix.
    /// \param lMax Values of l GREATER THAN this are included in the fiducial matrix.
    /// \param fwhm The full width at half maximum of the gaussian beam.
    /// \param goodPixels A pointer to a vector containing the indices of unmasked pixels, NULL to use them all.
    /// \param lp A pointer to a container of Legendre polynomials, NULL makes it calculate them by itself.
    /// \return A pointer to the generated fiducial matrix. The units are mK. It must be deleted after using.
    static CMatrix* getFiducialMatrix(const char* clFileName, long nSide, int lMax, double fwhm, const std::vector<int>* goodPixels = NULL, const LegendrePolynomialContainer* lp = NULL);
    
    /// Creates a covariance matrix corresponding to white noise.
    
    /// This function generates a diagonal matrix with fixed white noise in each pixel.
    /// \param nSide NSide of the matrix generated.
    /// \param noise The white noise in units of mK.
    /// \return A pointer to the generated noise matrix. The units are mK. It must be deleted after using.
    static CMatrix* generateNoiseMatrix(long nSide, double noise = 1e-3);
    
    /// Downgrade noise matrix from higher resolution to low resolution. Needs further testing!
    
    /// This function estimates the low resolution noise matrix using simulations. In high resolution the noise per pixel is determined from N_Obs. 
    /// Multiple instances of noise are generated in high resolution, then downgraded and multiplied by the beam.
    /// The noise matrix is estimated by calculating the correlations between different pixels in simulations.
    /// \param maskFileName The file containing the mask. The generated noise matrix is masked.
    /// \param noiseDataFileName This file contains the original maps in text format. Each line ust have pixel number, temperature, N_Obs separated by white space.
    /// \param sigma0 Sigma_0 of the original map in units of mK.
    /// \param fwhm Full width at half maximum of the beam for the low resolution map.
    /// \param nSideOriginal NSide of the high resolution map (NSide of low resolution is determined from the mask).
    /// \param fwhmOriginal Full width at half maximum of the beam of the original map.
    /// \return A pointer to the calculated noise matrix. The units are mK. It must be deleted after using.
    static CMatrix* calculateNoiseMatrix(const char* maskFileName, const char* noiseDataFileName, double sigma0, double fwhm, long nSideOriginal = 512, double fwhmOriginal = 1);
};

#endif
