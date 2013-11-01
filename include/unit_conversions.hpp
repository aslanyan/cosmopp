#ifndef COSMO_CPP_UNIT_CONVERSIONS_HPP
#define COSMO_CPP_UNIT_CONVERSIONS_HPP

#include <phys_constants.hpp>

namespace Phys
{
	
/// Convert Megaparsec to m.
/// \param x Distance in Mpc.
/// \return Distance in m.
inline
double MpcToM(double x)
{
	return x * MegaParsec;
}
	
/// Convert ev to Kg.
/// \param x Energy in ev.
/// \return Energy (or rather mass) in kg.
inline
double evToKg(double x)
{
	return x * eCharge / (cLight * cLight);
}
	
/// Convert inverse second to ev.
/// \param x Energy in inverse seconds.
/// \return Energy in ev.
inline
double inverseSecToEv(double x)
{
	return x * hbar / eCharge;
}
	
/// Convert inverse second to kg.
/// \param x Energy in inverse seconds.
/// \return Energy (or rather mass) in kg.
inline
double inverseSecToKg(double x)
{
	return x * hbar / (cLight * cLight);
}
	
/// Convert seconds to years.
/// \param x Time in seconds.
/// \return Time in years.
inline
double secToYear(double x)
{
	return x / 31556926;
}
	
/// Convert unitless to seconds.
/// Unitless quantities are in units hbar = c = kB = mPl = 1.
/// \param x Unitless quantity.
/// \return The quantity in seconds.
inline
double unitlessToSec(double x)
{
	return x * hbar / (mPl * cLight * cLight);
}
	
/// Convert seconds to unitless.
/// Unitless quantities are in units hbar = c = kB = mPl = 1.
/// \param x Quantity in seconds.
/// \return The quantity without units.
inline
double secToUnitless(double x)
{
	return x * (mPl * cLight * cLight) / hbar;
}

/// Convert unitless to inverse seconds.
/// Unitless quantities are in units hbar = c = kB = mPl = 1.
/// \param x Unitless quantity.
/// \return The quantity in inverse seconds.
inline
double unitlessToInverseSec(double x)
{
    return x * (mPl * cLight * cLight) / hbar;
}

/// Convert inverse seconds to unitless.
/// Unitless quantities are in units hbar = c = kB = mPl = 1.
/// \param x Quantity in inverse seconds.
/// \return The quantity without units.
inline
double inverseSecToUnitless(double x)
{
    return x * hbar / (mPl * cLight * cLight);
}
    
/// Convert meters to unitless.
/// Unitless quantities are in units hbar = c = kB = mPl = 1.
/// \param x Distance in meters.
/// \return Unitless distance.
inline
double mToUnitless(double x)
{
    return secToUnitless(x / cLight);
}
    
/// Convert unitless to meters.
/// Unitless quantities are in units hbar = c = kB = mPl = 1.
/// \param x Unitless distance.
/// \return Distance in meters.
inline
double unitlessToM(double x)
{
    return unitlessToSec(x) * cLight;
}

/// Convert Kelvin to unitless.
/// Unitless quantities are in units hbar = c = kB = mPl = 1.
/// \param x Temperature in K.
/// \return Unitless temperature.
inline
double kelvinToUnitless(double x)
{
    return kB * x / (mPl * cLight * cLight);
}
	
} //namespace Phys

#endif
