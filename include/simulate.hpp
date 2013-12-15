#ifndef COSMO_PP_SIMULATE_HPP
#define COSMO_PP_SIMULATE_HPP

#include <ctime>
#include <vector>

#include <whole_matrix.hpp>

#include <gmd.h>
#include <lavd.h>

#include <alm.h>
#include <xcomplex.h>
#include <healpix_map.h>

/// Simulation functions.

/// The functions in this class are used for simulating maps
class Simulate
{
public:
    /// Simulates alm from a whole matrix.
    
    /// This function generates a random instance of alm-s for which the covariance matrix is the given WholeMatrix.
    /// The generated alm-s correspond to a real map, so only m >= 0 are generated.
    /// \param wholeMatrix The covariance matrix in l-m space.
    /// \param alm The simulated alm to be returned.
    /// \param chi2 A pointer to chi2 that is to be returned. This is used for testing purposes only. If NULL then chi2 is not calculated.
    /// \param dof The number of degrees of freedom to be returned. Should be not NULL only if chi2 is not NULL. This is used for testing too.
    /// \param seed A random seed. If 0 then the current time (in seconds) is taken as the seed.
    static void simulateAlm(const WholeMatrix& wholeMatrix, Alm<xcomplex<double> >& alm, double* chi2 = NULL, int* dof = NULL, time_t seed = 0);

    /// Simulates alm from C_l-s.
    
    /// This function generates a random instance of alm-s with given C_l-s.
    /// The generated alm-s correspond to a real map, so only m >= 0 are generated.
    /// \param cl A vector of C_l-s, the index is l.
    /// \param alm The simulated alm to be returned.
    /// \param lMax The maximum value of l to use for simulation. In the default case of 0 it is determined from the size of cl.
    /// \param seed A random seed. If 0 then the current time (in seconds) is taken as the seed.
    static void simulateAlm(const std::vector<double>& cl, Alm<xcomplex<double> >& alm, int lMax = 0, time_t seed = 0);

    /// Diagonalize a square positive definite matrix.
    
    /// This function finds the eigenvalues of a matrix, then checks (if CHECKS_ON is defined) that the eigenvalues are real and positive, and the eigenvectors are orthonormal.
    /// \param matrix The matrix to be diagonalized.
    /// \param eigenvalsRe The vector to write the real parts of eigenvalues in.
    /// \param eigenvalsIm The vector to write the imaginary parts of eigenvalues in.
    /// \param vecs The matrix to write the eigenvectors in.
    static void diagonalizeMatrix(const LaGenMatDouble& matrix, LaVectorDouble& eigenvalsRe, LaVectorDouble& eigenvalsIm, LaGenMatDouble& vecs);

    /// White noise map generator.
    
    /// This function simulates a white noise map.
    /// \param map The map in which the simulated white noise will be written. The map's n_side or ordering scheme is not changed.
    /// \param noiseVal The value of the noise (default = 1.0).
    /// \param seed A random seed. If 0 then the current time (in seconds) is taken as the seed.
    static void simulateWhiteNoise(Healpix_Map<double>& map, double noiseVal = 1.0, time_t seed = 0);
};

#endif
