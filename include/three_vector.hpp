#ifndef COSMO_PP_THREE_VECTOR_HPP
#define COSMO_PP_THREE_VECTOR_HPP

#include <cmath>

#include <angular_coordinates.hpp>
#include <math_constants.hpp>

namespace Math
{

/// A three dimensional vector class.
template<typename T>
class ThreeVector
{
public:
    /// The variable type.
	typedef T VariableType;
	
public:
    /// Constructor.
    /// \param x The x coordinate.
    /// \param y The y coordinate.
    /// \param z The z coordinate.
	ThreeVector(VariableType x = 0, VariableType y = 0, VariableType z = 0) : x_(x), y_(y), z_(z) {}

    /// Destructor.
	~ThreeVector(){}
	
    /// Calculate norm squared.
    /// \return The norm squared.
	VariableType normSquared() const {return x_ * x_ + y_ * y_ + z_ * z_;}

    /// Dot product.
    /// \param other The other vector to take the dot product with.
    /// \return The dot product.
    VariableType operator * (const ThreeVector<T>& other) const {return x_ * other.x_ + y_ * other.y_ + z_ * other.z_;}
    
    /// Get the angular coordinates of the vector.
    /// \return The angular coordinates.
	AngularCoordinates angularCoordinates() const;
    
    /// Access x.
    /// \return The x value.
	VariableType x() const {return x_;}

    /// Access y.
    /// \return The y value.
	VariableType y() const {return y_;}

    /// Access z.
    /// \return The z value.
	VariableType z() const {return z_;}
    
    /// Access x.
    /// \return A reference to the x value (can be changed).
	VariableType& x() {return x_;}
    
    /// Access y.
    /// \return A reference to the y value (can be changed).
	VariableType& y() {return y_;}
    
    /// Access z.
    /// \return A reference to the z value (can be changed).
	VariableType& z() {return z_;}
    
    /// Set the coordinates.
    /// \param x The new x value.
    /// \param y The new y value.
    /// \param z The new z value.
    void set(VariableType x, VariableType y, VariableType z) { x_ = x; y_ = y; z_ = z; }
	
    /// Copy form another vector.
    /// \param v The vector to copy from.
    /// \return A reference to this after copying.
	ThreeVector<T>& operator = (const ThreeVector<T>& v) { x_ = v.x_; y_ = v.y_; z_ = v.z_; return *this; }
	
private:
	VariableType x_;
	VariableType y_;
	VariableType z_;
};

typedef ThreeVector<double> ThreeVectorDouble;

template<typename T>
AngularCoordinates ThreeVector<T>::angularCoordinates() const
{
	AngularCoordinates res;
	double norm = std::sqrt(double(normSquared()));
	if(norm == 0)
		return res;
	
	res.theta = std::acos(z_ / norm);
	norm = std::sqrt(double(x_ * x_ + y_ * y_));
	if(norm == 0)
		return res;
	
	res.phi = std::acos(x_ / norm);
	if(y_ < 0)
		res.phi = 2 * pi - res.phi;
	
	return res;
}
	
} //namespace Math

#endif
