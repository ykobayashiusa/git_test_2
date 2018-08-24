/*
 * Copyright (C) 2017 Matthias Kirchhart
 *
 * This file is part of vorticus.
 *
 * vorticus is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3, or (at your option) any later
 * version.
 *
 * vorticus is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * vorticus; see the file COPYING.  If not see http://www.gnu.org/licenses.
 */

namespace geometry
{

///////////////////////
// Member functions. //
///////////////////////

/*!
 * \brief Construct from given cartesian coordinates.
 */
inline
tl_point::tl_point( point p ) noexcept:
 x { p.x }, y { p.y }, z { p.z }
{}

/*!
 * \brief Construct from given cartesian coordinates.
 */
inline
tl_point::tl_point( tl_real xx, tl_real yy, tl_real zz ) noexcept:
 x { xx }, y { yy }, z { zz }
{}

/*!
 * \brief Construct from given cartesian coordinates.
 */
inline
tl_point::tl_point( point p, tensor dx ) noexcept:
 x { p.x, std::array<real,3> {{ dx(0,0), dx(0,1), dx(0,2) }} }, 
 y { p.y, std::array<real,3> {{ dx(1,0), dx(1,1), dx(1,2) }} }, 
 z { p.z, std::array<real,3> {{ dx(2,0), dx(2,1), dx(2,2) }} }
{}


/*!
 * \brief Returns the Euclidean length of the point.
 */
inline
tl_real tl_point::r() const noexcept
{
	return sqrt( x*x + y*y + z*z );
}

/*!
 * \brief Returns the squared Euclidean length of the point.
 */
inline
tl_real tl_point::r2() const noexcept
{
	return x*x + y*y + z*z;
}

/*!
 * \brief Unary minus.
 */
inline
tl_point tl_point::operator-() const noexcept
{
	return tl_point { -x, -y, -z };
}

/*!
 * \brief Vector addition.
 */
inline
tl_point tl_point::operator+( const point& rhs ) const noexcept
{
	return tl_point { x + rhs.x, y + rhs.y, z + rhs.z };
}

/*!
 * \brief Vector substraction.
 */
inline
tl_point tl_point::operator-( const point& rhs ) const noexcept
{
	return tl_point { x - rhs.x, y - rhs.y, z - rhs.z };
}

/*!
 * \brief Vector addition.
 */
inline
tl_point tl_point::operator+( const tl_point& rhs ) const noexcept
{
	return tl_point { x + rhs.x, y + rhs.y, z + rhs.z };
}

/*!
 * \brief Vector substraction.
 */
inline
tl_point tl_point::operator-( const tl_point& rhs ) const noexcept
{
	return tl_point { x - rhs.x, y - rhs.y, z - rhs.z };
}

/*!
 * \brief Vector addition.
 */
inline
void tl_point::operator+=( const point& rhs ) noexcept
{
	x += rhs.x;
	y += rhs.y;
	z += rhs.z;
}

/*!
 * \brief Vector substraction.
 */
inline
void tl_point::operator-=( const point& rhs ) noexcept
{
	x -= rhs.x;
	y -= rhs.y;
	z -= rhs.z;
}

/*!
 * \brief Vector addition.
 */
inline
void tl_point::operator+=( const tl_point& rhs ) noexcept
{
	x += rhs.x;
	y += rhs.y;
	z += rhs.z;
}

/*!
 * \brief Vector substraction.
 */
inline
void tl_point::operator-=( const tl_point& rhs ) noexcept
{
	x -= rhs.x;
	y -= rhs.y;
	z -= rhs.z;
}

/*!
 * \brief Multiplication with a scalar.
 */
inline
tl_point tl_point::operator*( real rhs ) const noexcept
{
	return tl_point { x*rhs, y*rhs, z*rhs };
}

/*!
 * \brief Division by a scalar.
 */
inline
tl_point tl_point::operator/( real rhs ) const noexcept
{
	return tl_point { x/rhs, y/rhs, z/rhs };
}

/*!
 * \brief Multiplication with a scalar.
 */
inline
tl_point tl_point::operator*( tl_real rhs ) const noexcept
{
	return tl_point { x*rhs, y*rhs, z*rhs };
}

/*!
 * \brief Division by a scalar.
 */
inline
tl_point tl_point::operator/( tl_real rhs ) const noexcept
{
	return tl_point { x/rhs, y/rhs, z/rhs };
}

/*!
 * \brief Multiplication with a scalar.
 */
inline
void tl_point::operator*=( real rhs ) noexcept
{
	x *= rhs;
	y *= rhs;
	z *= rhs;
}

/*!
 * \brief Division by a scalar.
 */
inline
void tl_point::operator/=( real rhs ) noexcept
{
	x /= rhs;
	y /= rhs;
	z /= rhs;
}

/*!
 * \brief Multiplication with a scalar.
 */
inline
void tl_point::operator*=( tl_real rhs ) noexcept
{
	x *= rhs;
	y *= rhs;
	z *= rhs;
}

/*!
 * \brief Division by a scalar.
 */
inline
void tl_point::operator/=( tl_real rhs ) noexcept
{
	x /= rhs;
	y /= rhs;
	z /= rhs;
}

inline
point tl_point::val() const noexcept
{
    return point { x.val, y.val, z.val };
}

inline
tensor tl_point::grad() const noexcept
{
    return tensor
    {
        x.derivatives[0], x.derivatives[1], x.derivatives[2],
        y.derivatives[0], y.derivatives[1], y.derivatives[2],
        z.derivatives[0], z.derivatives[1], z.derivatives[2]
    };
}

inline
point tl_point::curl() const noexcept
{
    return point
    {
        z.derivatives[1] - y.derivatives[2],
        x.derivatives[2] - z.derivatives[0],
        y.derivatives[0] - x.derivatives[1]
    };
}

///////////////////////////
// Non-member functions. //
///////////////////////////

/*!
 * \brief Returns the length of the given point.
 */
inline
tl_real length( const tl_point& p ) noexcept
{
	return p.r();
}

/*!
 * \brief Multiplication by a scalar.
 */
inline
tl_point operator+( const point& lhs, const tl_point &rhs ) noexcept
{
	return rhs + lhs;
}

/*!
 * \brief Multiplication by a scalar.
 */
inline
tl_point operator-( const point& lhs, const tl_point &rhs ) noexcept
{
	return -(rhs - lhs);
}


/*!
 * \brief Multiplication by a scalar.
 */
inline
tl_point operator*( const tl_real &lhs, const tl_point& rhs ) noexcept
{
	return rhs*lhs;
}

/*!
 * \brief Multiplication by a scalar.
 */
inline
tl_point operator*( const real lhs, const tl_point& rhs ) noexcept
{
	return rhs*lhs;
}

/*!
 * \brief Multiplication by a scalar.
 */
inline
tl_point operator*( const tl_real &lhs, const point rhs ) noexcept
{
	return tl_point { lhs*rhs.x, lhs*rhs.y, lhs*rhs.z };
}

/*!
 * \brief Returns the scalar product of the two given points.
 */
inline
tl_real scal_prod( const tl_point& lhs, const tl_point& rhs ) noexcept
{
	return lhs.x*rhs.x +
	       lhs.y*rhs.y +
	       lhs.z*rhs.z;
}

/*!
 * \brief Returns the scalar product of the two given points.
 */
inline
tl_real scal_prod( const tl_point& lhs, const point& rhs ) noexcept
{
	return lhs.x*rhs.x +
	       lhs.y*rhs.y +
	       lhs.z*rhs.z;
}

/*!
 * \brief Returns the scalar product of the two given points.
 */
inline
tl_real scal_prod( const point& lhs, const tl_point& rhs ) noexcept
{
	return rhs.x*lhs.x +
	       rhs.y*lhs.y +
	       rhs.z*lhs.z;
}

/*!
 * \brief Returns the cross product of the given points.
 */
inline
tl_point cross_prod( const tl_point& lhs, const tl_point& rhs ) noexcept
{
	return tl_point{ lhs.y*rhs.z - lhs.z*rhs.y,
	                 lhs.z*rhs.x - lhs.x*rhs.z,
	                 lhs.x*rhs.y - lhs.y*rhs.x };
}

/*!
 * \brief Returns the cross product of the given points.
 */
inline
tl_point cross_prod( const point& lhs, const tl_point& rhs ) noexcept {
	return  tl_point{ rhs.z*lhs.y - rhs.y*lhs.z,
	                  rhs.x*lhs.z - rhs.z*lhs.x,
	                  rhs.y*lhs.x - rhs.x*lhs.y };
}

/*!
 * \brief Returns the cross product of the given points.
 */
inline
tl_point cross_prod( const tl_point& lhs, const point& rhs ) noexcept
{
	return tl_point{ lhs.y*rhs.z - lhs.z*rhs.y,
	                 lhs.z*rhs.x - lhs.x*rhs.z,
	                 lhs.x*rhs.y - lhs.y*rhs.x };
}

}

