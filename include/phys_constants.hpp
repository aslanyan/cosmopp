#ifndef COSMO_CPP_PHYS_CONSTANTS_HPP
#define COSMO_CPP_PHYS_CONSTANTS_HPP

#include <math_constants.hpp>

/// @namespace Phys
/// Contains general physics related functions, classes, and constants.
namespace Phys
{
	
/// Speed of light (m / s)
const double cLight = 299792458;

/// Planck's constant (J s)
const double hPlanck = 6.62606896E-34;

/// Barred Planck's constant (J s)
const double hbar = hPlanck / (2 * Math::pi);
	
/// Gravitational constant (m^3 kg^{-1} s^{-2})
const double GNewton = 6.67428E-11;
	
/// Reduced Planck mass (kg)
const double mPl = 4.34136451390282E-9;
	
/// Elementary charge (C)
const double eCharge = 1.602176487E-19;

/// Boltzmann's constant (m^2 kg s^{-2} K^{-1})
const double kB = 1.3806504E-23;
	
/// 1 Megaparsec (m)
const double MegaParsec = 3.085677581282e22;

} //namespace Phys

#endif
