#ifndef COSMO_PP_INTEGRAL_HPP
#define COSMO_PP_INTEGRAL_HPP

#include <algorithm>

#include <function.hpp>
#include <macros.hpp>

/// @namespace Math
/// Contains math related functions, classes, and constants.
namespace Math
{

/// One dimensional integral.

/// Numerically calculates the 1-dimensional integral of a function using the trapezoids method on equal length intervals.
/// \param f The function to integrate.
/// \param xMin The lower limit of the interval.
/// \param xMax The upper limit of the interval.
/// \param numberOfPoints The number of points to split the interval into.
/// \return The value of the integral.
inline
double realIntegral1D(const RealFunction& f, double xMin, double xMax, int numberOfPoints = 1000)
{
	int sign = 1;
	if(xMax == xMin)
		return 0;
	
	if(xMax < xMin)
	{
		std::swap(xMin, xMax);
		sign = -1;
	}
	
	check(numberOfPoints > 1, "");
	
	const double delta = (xMax - xMin) / (numberOfPoints - 1);
	double a = xMin, b = xMin + delta;
	double f1 = f.evaluate(a), f2 = f.evaluate(b);
	double r = (f1 + f2) / 2;
	for(int i = 1; i < numberOfPoints - 2; ++i)
	{
		a = b;
		b = a + delta;
		f1 = f2;
		f2 = f.evaluate(b);
		r += (f1 + f2) / 2;
	}
	a = b;
	b = xMax;
	f1 = f2;
	f2 = f.evaluate(b);
	r += (f1 + f2) / 2;
	
	return r * delta * sign;
}
	
} //namespace Math

#endif
