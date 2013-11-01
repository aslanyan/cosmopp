#ifndef COSMO_CPP_INT_OPERATIONS_HPP
#define COSMO_CPP_INT_OPERATIONS_HPP

namespace Math
{

/// Integer square root function.

/// A fast way to take the square root of an integer.
/// \param n The integer to take the square root of.
/// \return The integer part of the square root of n.
inline
unsigned int isqrt(unsigned long n)
{
	unsigned int c = 0x8000;
	unsigned int g = 0x8000;
	
	for(;;) {
		if(g*g > n)
			g ^= c;
		c >>= 1;
		if(c == 0)
			return g;
		g |= c;
	}
}

} //namespace Math

#endif

