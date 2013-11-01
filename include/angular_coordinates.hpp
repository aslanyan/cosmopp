#ifndef COSMO_CPP_ANGULAR_COORDINATES_HPP
#define COSMO_CPP_ANGULAR_COORDINATES_HPP

namespace Math
{

/// Angular Coordinates.
struct AngularCoordinates
{
    /// Constructor.
    /// \param t Angular coordinate theta
    /// \param p Angular coordinate phi
	AngularCoordinates(double t = 0, double p = 0) : theta(t), phi(p) {}
	
    /// Angular coordinate theta
	double theta;

    /// Angular coordinate phi
	double phi;
};
	
} //namespace Math

#endif
