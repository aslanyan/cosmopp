#ifndef COSMO_CPP_THREE_ROTATION_HPP
#define COSMO_CPP_THREE_ROTATION_HPP

#include <cmath>
#include <vector>

#include <macros.hpp>
#include <three_vector.hpp>
#include <numerics.hpp>

namespace Math
{

/// A class for rotation matrix in three dimensions, i.e. an element of O(3).

/// This Matrix as defined does a PASSIVE (i.e. coordinate transformation) as follows:
/// First rotate counterclockwise around the z axis by angle phi
/// then rotate counterclockwise around the new x axis by angle theta
/// then rotate counterclockwise around the new z axis by angle psi
/// Thus we obtain a new coordinate frame, applying the rotation to a three-vector will give
/// the same vector in the new frame.
class ThreeRotationMatrix
{
public:
	typedef std::vector<double> VectorType;
	
public:
    /// Constructor.
    /// \param phi Euler angle phi.
    /// \param theta Euler angle theta.
    /// \param psi Euler angle psi.
	ThreeRotationMatrix(double phi = 0, double theta = 0, double psi = 0);

    /// Copy constructor.
    /// \param other The matrix to copy from.
	ThreeRotationMatrix(const ThreeRotationMatrix& other);
	
    /// Destructor.
	~ThreeRotationMatrix() {}
	
    /// Set the Euler angles to new values.
    /// \param phi Euler angle phi.
    /// \param theta Euler angle theta.
    /// \param psi Euler angle psi.
	void set(double phi, double theta, double psi);
	
    /// Set the rotation matrix to be around a given axis by a given angle.
    /// \param axis The direction of the axis.
    /// \param phi The rotation angle.
	void set(ThreeVectorDouble axis, double phi);
	
    /// Access to a given row.
    /// \param i The index of the row (between 0 and 2).
    /// \return A constant reference to the vector representing row i.
	const VectorType& operator[](int i) const { check(i >= 0 && i < 3, "invalid index"); return matrix_[i]; }

    /// Apply the matrix to a given three-vector.
    /// \param v The vector to apply the rotation to.
    /// \return The rotated three-vector.
	ThreeVectorDouble operator * (const ThreeVectorDouble& v) const;
	
    /// Multiply by another rotation matrix. 
    /// The order of multiplication is this * other.
    /// \param other The other rotation matrix.
    /// \return The resulting rotation matrix after the multiplication. Note that this isn't changed.
	ThreeRotationMatrix operator * (const ThreeRotationMatrix& other) const;

    /// Copy operator.
    /// \param other The rotation matrix to copy from.
    /// \return A reference to this after the copy is done.
	ThreeRotationMatrix& operator = (const ThreeRotationMatrix& other);
	
    /// Check equality.
    /// The precision of equality is 1e-3.
    /// \param other The other matrix to compare to.
    /// \return The result of the comparison.
	bool operator == (const ThreeRotationMatrix& other) const;
	
private:
	std::vector<VectorType> matrix_;
};
	
inline
ThreeRotationMatrix::ThreeRotationMatrix(double phi, double theta, double psi)
{
	matrix_.resize(3);
	for(int i = 0; i < 3; ++i)
		matrix_[i].resize(3);
	
	set(phi, theta, psi);
}
	
inline
ThreeRotationMatrix::ThreeRotationMatrix(const ThreeRotationMatrix& other)
{
	ThreeRotationMatrix();
	*this = other;
}

inline
void ThreeRotationMatrix::set(double phi, double theta, double psi)
{
	const double c1 = std::cos(phi), s1 = std::sin(phi);
	const double c2 = std::cos(theta), s2 = std::sin(theta);
	const double c3 = std::cos(psi), s3 = std::sin(psi);
	
	matrix_[0][0] = c1 * c3 - c2 * s1 * s3;
	matrix_[0][1] = c2 * c1 * s3 + c3 * s1;
	matrix_[0][2] = s3 * s2;
	matrix_[1][0] = -c1 * s3 - c3 * c2 * s1;
	matrix_[1][1] = c1 * c2 * c3 - s1 * s3;
	matrix_[1][2] = c3 * s2;
	matrix_[2][0] = s2 * s1;
	matrix_[2][1] = -c1 * s2;
	matrix_[2][2] = c2;
}
	
inline
void ThreeRotationMatrix::set(ThreeVectorDouble axis, double phi)
{
	const double normSquared = axis.normSquared();
	check(normSquared != 0, "axis must be nonzero");
	const double c = std::cos(phi);
	const double s = std::sin(phi);
	const double& u = axis.x();
	const double& v = axis.y();
	const double& w = axis.z();
	const double l = std::sqrt(normSquared);
	
	matrix_[0][0] = (u * u + (v * v + w * w) * c) / normSquared;
	matrix_[0][1] = (u * v * (1 - c) + w * l * s) / normSquared;
	matrix_[0][2] = (u * w * (1 - c) - v * l * s) / normSquared;
	matrix_[1][0] = (u * v * (1 - c) - w * l * s) / normSquared;
	matrix_[1][1] = (v * v + (u * u + w * w) * c) / normSquared;
	matrix_[1][2] = (v * w * (1 - c) + u * l * s) / normSquared;
	matrix_[2][0] = (u * w * (1 - c) + v * l * s) / normSquared;
	matrix_[2][1] = (v * w * (1 - c) - u * l * s) / normSquared;
	matrix_[2][2] = (w * w + (u * u + v * v) * c) / normSquared;
}

/*
inline
void ThreeRotationMatrix::setRotationAroundX(double phi)
{
	set(0, phi, 0);
}
	
inline
void ThreeRotationMatrix::setRotationAroundY(double phi)
{
	matrix_[0][0] = std::cos(phi);
	matrix_[0][1] = 0;
	matrix_[0][2] = std::sin(phi);
	matrix_[1][0] = 0;
	matrix_[1][1] = 1;
	matrix_[1][2] = 0;
	matrix_[2][0] = -matrix_[0][2];
	matrix_[2][1] = 0;
	matrix_[2][2] = matrix_[0][0];
}
	
inline
void ThreeRotationMatrix::setRotationAroundZ(double phi)
{
	set(phi, 0, 0);
}
*/
	
inline
ThreeVectorDouble ThreeRotationMatrix::operator * (const ThreeVectorDouble& v) const
{
	const double x = matrix_[0][0] * v.x() + matrix_[0][1] * v.y() + matrix_[0][2] * v.z();
	const double y = matrix_[1][0] * v.x() + matrix_[1][1] * v.y() + matrix_[1][2] * v.z();
	const double z = matrix_[2][0] * v.x() + matrix_[2][1] * v.y() + matrix_[2][2] * v.z();
	
	const ThreeVectorDouble res(x, y, z);
	check((res.normSquared() == 0 && v.normSquared() == 0) || std::abs(res.normSquared() - v.normSquared()) / v.normSquared() < 1E-5,
		  "rotation should preserve the norm");
	
	return res;
}
	
inline
ThreeRotationMatrix ThreeRotationMatrix::operator * (const ThreeRotationMatrix& other) const
{
	ThreeRotationMatrix result;
	check(matrix_.size() == 3, "");
	check(other.matrix_.size() == 3, "");
	check(result.matrix_.size() == 3, "");
	for(int i = 0; i < 3; ++i)
	{
		check(matrix_[i].size() == 3, "");
		check(other.matrix_[i].size() == 3, "");
		check(result.matrix_[i].size() == 3, "");
		
		for(int j = 0; j < 3; ++j)
		{
			result.matrix_[i][j] = 0;
			for(int k = 0; k < 3; ++k)
				result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
		}
	}
	
	return result;
}

inline
ThreeRotationMatrix& ThreeRotationMatrix::operator = (const ThreeRotationMatrix& other)
{
	check(matrix_.size() == 3, "");
	check(other.matrix_.size() == 3, "");
	for(int i = 0; i < 3; ++i)
	{
		check(matrix_[i].size() == 3, "");
		check(other.matrix_[i].size() == 3, "");
		for(int j = 0; j < 3; ++j)
			matrix_[i][j] = other.matrix_[i][j];
	}
	return *this;
}
	
inline
bool ThreeRotationMatrix::operator == (const ThreeRotationMatrix& other) const
{
	check(matrix_.size() == 3, "");
	check(other.matrix_.size() == 3, "");
	for(int i = 0; i < 3; ++i)
	{
		check(matrix_[i].size() == 3, "");
		check(other.matrix_[i].size() == 3, "");
		for(int j = 0; j < 3; ++j)
		{
			if(!areEqual(matrix_[i][j], other.matrix_[i][j], 1E-3))
				return false;
		}
	}
	return true;
}
	
} //namespace Math

#endif
