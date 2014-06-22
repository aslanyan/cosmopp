#ifndef COSMO_PP_UTILS_HPP
#define COSMO_PP_UTILS_HPP

#include <vector>

#include <healpix_map.h>

/// A class for various utility functions.

/// All of the functions in the class are static.
class Utils
{
public:
    /// Reads a mask file.

    /// This function reads a mask file. A mask is a healpix map with pixels equal to 0 or 1.
    /// \param maskFileName The file of a fits file containing the mask.
    /// \param nSide NSide of the mask is written here.
    /// \param goodPixels A vector that upon return contains the original indices of the unmasked pixels (i.e. 1).
    static void readMask(const char* maskFileName, long& nSide, std::vector<int>& goodPixels);
    
    /// Calculates the gaussian beam function.
    
    /// This function calculates the gaussian beam function coefficien B_l
    /// \param l The l value.
    /// \param fwhm Full width at half maximum of the gaussian beam in degrees.
    /// \return B_l.
    static double beamFunction(int l, double fwhm);
    
    /// Reads a pixel window function.
    
    /// Reads a pixel window function from the Healpix data directory, also multiplies by a gaussian beam function.
    /// \param f The gaussian beam coefficiants vector to be returned, the index is l.
    /// \param nSide The NSide of the pixel window function.
    /// \param lMax Maximum l value to be read.
    /// \param fwhm The full width at half maximum of the gaussian beam in degrees. 0 means no beam.
    /// \param polarization Specifies if polarization pixel window function should be read instead of temperature (false by default).
    static void readPixelWindowFunction(std::vector<double>& f, long nSide, int lMax, double fwhm = 0, bool polarization = false);

    /// Reads the values of cl-s from a text file.

    /// Reads the values of cl-s from a text file. Each line should contain one value, starting from l = 0. The l values can proceed the cl values, in this case this should be indicated by the parameter hasL.
    /// \param fileName The name of the file.
    /// \param cl The vector where the cl values will be written.
    /// \param hasL Indicates if the l values proceed the cl values, false by default.
    /// \param isDl Indicates if the cl-s include a factor of l(l+1)/2pi. If so, this will be corrected for in the output (false by default).
    static void readClFromFile(const char* fileName, std::vector<double>& cl, bool hasL = false, bool isDl = false);

    /// Mask out a circular region around a given direction.
    /// \param mask The mask to be modified.
    /// \param theta The theta angle of the direction of the center of the region.
    /// \param phi The phi angle of the direction of the center of the region.
    /// \param angle The radius of the circular region.
    static void maskRegion(Healpix_Map<double>& mask, double theta, double phi, double angle);

    /// Mask out multiple circular regions around given directions.
    /// \param mask The mask to be modified.
    /// \param theta A vector of the theta angles of the directions of the centers of the regions.
    /// \param phi A vector of the phi angles of the directions of the centers of the regions. Must have the same size as theta.
    /// \param angle A vector of the radii of the circular regions. Must have the same size as theta.
    static void maskRegions(Healpix_Map<double>& mask, const std::vector<double>& theta, const std::vector<double>& phi, const std::vector<double>& angle);
};

#endif
