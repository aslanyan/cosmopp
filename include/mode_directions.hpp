#ifndef COSMO_CPP_MODE_DIRECTIONS_HPP
#define COSMO_CPP_MODE_DIRECTIONS_HPP

#include <alm.h>
#include <xcomplex.h>
#include <healpix_map.h>

/// A class for calculating angular momentum dispersion.
class ModeDirections
{
public:
    /// Constructor.
    /// \param alm The alm-s of the map under consideration.
    ModeDirections(const Alm<xcomplex<double> >& alm);
    
    /// Calculate the angular momentum dispersion for a given mode in a given direction.
    /// \param l The mode l.
    /// \param theta The theta of the direction.
    /// \param phi The phi of the direction.
    /// \return The angular momentum dispersion in the given direction.
    double calculateAngularMomentumDispersion(int l, double theta, double phi) const;

    /// Find the direction in which the angular momentum dispersion is maximized. The maximization is done by iterating over the sphere with Healpix pixelization.
    /// \param l The mode l.
    /// \param nSide The nSide which determines the pixelization to use for finding the maximum. Higher values will give more accurate results but will slow down the calculation.
    /// \param theta The angle theta of the direction maximizing the angular momentum dispersion will be written here.
    /// \param phi The angle phi of the direction maximizing the angular momentum dispersion will be written here.
    /// \param map A pointer to a Healpix map where the angular momentum dispersion as a function of direction will be stored. Give NULL to not store it (this is the default option). The resulting N_side of the map is the same as nSide.
    void maximizeAngularMomentumDispersion(int l, long nSide, double& theta, double& phi, Healpix_Map<double>* map = NULL) const;
    
private:
    Alm<xcomplex<double> > alm_;
};

#endif
