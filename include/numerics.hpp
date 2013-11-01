#ifndef COSMO_CPP_NUMERICS_HPP
#define COSMO_CPP_NUMERICS_HPP

#include <cmath>

namespace Math
{

/// Checks if two numbers are equal with a given precision.

/// The real type numbers that should be equal could differ slightly because of numerics. This function checks their equality with a given precision.
/// \param a The first number.
/// \param b The second number.
/// \param precision The precision with which to check the numbers equality.
/// \return true if b is equal to a with the given precision.
template<typename T>
bool areEqual(T a, T b, T precision = 1E-7)
{
	if(std::abs(a) < precision)
		return std::abs(b) < precision;
	return std::abs(a - b) / std::abs(a) < precision;
}
	
/// This class is used to check if a real number is less than the other one with a given precision. 
/// It can be used to construct sets and maps of real numbers. If the numbers are equal within a given precision then they're not considered to be less than one another.
template<typename T>
class Less
{
public:
    /// Constructor.
    /// \param precision The precision within which the numbers should be considered equal.
	Less(T precision = 1E-7) : precision_(precision) {}
	
    /// The less operator.
    /// \param a First number.
    /// \param b Second number.
    /// \return true if a < b given the fact that they're considered equal if they're within the precision of each other.
	bool operator() (T a, T b) const
	{
		if(areEqual(a, b, precision_))
			return false;
		
		return a < b;
	}
	
private:
	T precision_;
};
	
} //namespace Math

#endif
