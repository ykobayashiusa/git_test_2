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

/*!
 * \file geometry/cmplx_point.tpp
 * \author Matthias Kirchhart <matthias@kirchhart.de>
 * \brief Implementation of the cmplx_point class.
 */

namespace geometry
{

///////////////////////
// Member functions. //
///////////////////////

constexpr
cmplx_point::cmplx_point( const point& rhs ) noexcept:
 x { rhs.x, 0 }, y { rhs.y, 0 }, z { rhs.z, 0 }
{}

constexpr
cmplx_point::cmplx_point( const point& real, const point& imag ) noexcept:
 x { real.x, imag.x }, y { real.y, imag.y }, z { real.z, imag.z }
{}

constexpr
cmplx_point::cmplx_point( cmplx xx, cmplx yy, cmplx zz ) noexcept:
 x { xx }, y { yy }, z { zz }
{}

inline
point cmplx_point::real() const noexcept
{
	return point { x.real(), y.real(), z.real() };
}

inline
point cmplx_point::imag() const noexcept
{
	return point { x.imag(), y.imag(), z.imag() };
}

inline
cmplx_point cmplx_point::conj() const noexcept
{
	return cmplx_point {  std::conj(x), std::conj(y), std::conj(z) };
}

inline
cmplx_point cmplx_point::operator+( const cmplx_point& rhs ) const noexcept
{
	return cmplx_point { x+rhs.x, y+rhs.y, z+rhs.z };
}

inline
cmplx_point cmplx_point::operator-( const cmplx_point& rhs ) const noexcept
{
	return cmplx_point { x-rhs.x, y-rhs.y, z-rhs.z };
}

inline
cmplx_point cmplx_point::operator-() const noexcept
{
	return cmplx_point { -x, -y, -z };
}

inline
void cmplx_point::operator+=( const cmplx_point& rhs ) noexcept
{
	x += rhs.x;
	y += rhs.y;
	z += rhs.z;
}

inline
void cmplx_point::operator-=( const cmplx_point& rhs ) noexcept
{
	x -= rhs.x;
	y -= rhs.y;
	z -= rhs.z;
}

inline
cmplx_point cmplx_point::operator*( cmplx rhs ) const noexcept
{
	return cmplx_point { x*rhs, y*rhs, z*rhs };
}

inline
cmplx_point cmplx_point::operator/( cmplx rhs ) const noexcept
{
	return cmplx_point { x/rhs, y/rhs, z/rhs };
}

inline
void cmplx_point::operator*=( cmplx rhs ) noexcept
{
	x *= rhs;
	y *= rhs;
	z *= rhs;
}

inline
void cmplx_point::operator/=( cmplx rhs ) noexcept
{
	x /= rhs;
	y /= rhs;
	z /= rhs;
}

///////////////////////////
// Non-member functions. //
///////////////////////////

inline
cmplx_point operator*( cmplx lhs, const cmplx_point& rhs ) noexcept
{
	return rhs*lhs;
}

inline
real length( const cmplx_point& p ) noexcept
{
	return std::sqrt(std::abs(scal_prod(p,p)));
}

inline
cmplx scal_prod( const cmplx_point& lhs, const cmplx_point& rhs ) noexcept
{
	return conj(lhs.x)*rhs.x +
	       conj(lhs.y)*rhs.y +
	       conj(lhs.z)*rhs.z;
}

inline
cmplx_point cross_prod( const cmplx_point& lhs, const cmplx_point& rhs ) noexcept
{
	return cmplx_point { lhs.y*rhs.z - lhs.z*rhs.y,
	                     lhs.z*rhs.x - lhs.x*rhs.z,
		                 lhs.x*rhs.y - lhs.y*rhs.x };
}

inline
cmplx_point conj( const cmplx_point& p ) noexcept
{
	return cmplx_point { conj(p.x), conj(p.y), conj(p.z) };
}

}

