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
 * \file geometry/point.tpp
 * \author Matthias Kirchhart <matthias@kirchhart.de>
 * \brief Implementation of the geometry::point class.
 */

namespace geometry
{

///////////////////////
// Member functions. //
///////////////////////

/*!
 * \brief Construct from given cartesian coordinates.
 */
constexpr
point::point( real xx, real yy, real zz ) noexcept:
 x { xx }, y { yy }, z { zz }
{}


/*!
 * \brief Returns the Euclidean length of the point.
 */
inline
real point::r() const noexcept
{
	return std::sqrt( x*x + y*y + z*z );
}

/*!
 * \brief Returns the squared Euclidean length of the point.
 */
constexpr
real point::r2() const noexcept
{
	return x*x + y*y + z*z;
}

/*!
 * \brief Returns the polar angle \f$\theta\f$ of the point.
 */
inline
real point::theta() const noexcept
{
	return std::acos( z / r() );
}

/*!
 * \brief Returns the azimuthal angle \f$\varphi\f$ of the point.
 */
inline
real point::phi() const noexcept
{
	return std::atan2( y, x );
}

/*!
 * \brief Returns the Manhattan-distance norm of the point.
 */
inline
real point::max_norm() const noexcept
{
    using std::max;
    using std::abs;
    return max( abs(x), max( abs(y), abs(z) ) );
}

/*!
 * \brief Returns first unit vector of spherical coordinate system.
 *
 * This method returns the unit vector which points to the direction of
 * increasing r.
 */
inline
point point::e_r() const noexcept
{
	return *this / r();
}

/*!
 * \brief Returns second unit vector of spherical coordinate system.
 *
 * This method returns the unit vector which points to the direction of
 * increasing \f$\theta\f$.
 */
inline
point point::e_theta() const noexcept
{
    using std::sin;
    using std::cos;
    return point { cos(theta())*cos(phi()),
                   cos(theta())*sin(phi()),
                   -sin(theta()) };
}

/*!
 * \brief Returns third unit vector of spherical coordinate system.
 *
 * This method returns the unit vector which points to the direction of
 * increasing \f$\varphi\f$.
 */
inline
point point::e_phi() const noexcept
{
    using std::sin;
    using std::cos;
	return point { -sin(phi()),
	                cos(phi()),
	                0 };
}

/*!
 * \brief Unary minus.
 */
constexpr
point point::operator-() const noexcept
{
	return point { -x, -y, -z };
}

/*!
 * \brief Vector addition.
 */
constexpr
point point::operator+( const point& rhs ) const noexcept
{
	return point { x + rhs.x, y + rhs.y, z + rhs.z };
}

/*!
 * \brief Vector substraction.
 */
constexpr
point point::operator-( const point& rhs ) const noexcept
{
	return point { x - rhs.x, y - rhs.y, z - rhs.z };
}

/*!
 * \brief Vector addition.
 */
inline
void point::operator+=( const point& rhs ) noexcept
{
	x += rhs.x;
	y += rhs.y;
	z += rhs.z;
}

/*!
 * \brief Vector substraction.
 */
inline
void point::operator-=( const point& rhs ) noexcept
{
	x -= rhs.x;
	y -= rhs.y;
	z -= rhs.z;
}

/*!
 * \brief Multiplication with a scalar.
 */
constexpr
point point::operator*( real rhs ) const noexcept
{
	return point { x*rhs, y*rhs, z*rhs };
}

/*!
 * \brief Division by a scalar.
 */
constexpr
point point::operator/( real rhs ) const noexcept
{
	return point { x/rhs, y/rhs, z/rhs };
}

/*!
 * \brief Multiplication with a scalar.
 */
inline
void point::operator*=( real rhs ) noexcept
{
	x *= rhs;
	y *= rhs;
	z *= rhs;
}

/*!
 * \brief Division by a scalar.
 */
inline
void point::operator/=( real rhs ) noexcept
{
	x /= rhs;
	y /= rhs;
	z /= rhs;
}

///////////////////////////
// Non-member functions. //
///////////////////////////

/*!
 * \brief Multiplication by a scalar.
 */
constexpr
point operator*( real lhs, const point& rhs ) noexcept
{
	return rhs*lhs;
}

/*!
 * \brief Returns the length of the given point.
 */
inline
real length( const point& p ) noexcept
{
	return p.r();
}

inline
point min_components( const point& lhs, const point& rhs ) noexcept
{
	return point { std::min( lhs.x, rhs.x ),
	               std::min( lhs.y, rhs.y ),
	               std::min( lhs.z, rhs.z ) };
}

inline
point max_components( const point& lhs, const point& rhs ) noexcept
{
	return point { std::max( lhs.x, rhs.x ),
	               std::max( lhs.y, rhs.y ),
	               std::max( lhs.z, rhs.z ) };
}

constexpr
bool operator==( const point& lhs, const point& rhs ) noexcept
{
    return lhs.x == rhs.x &&
           lhs.y == rhs.y &&
           lhs.z == rhs.z;
}

constexpr
bool operator!=( const point& lhs, const point& rhs ) noexcept
{
    return ! ( lhs == rhs );
}

/*!
 * \brief Returns the scalar product of the two given points.
 */
constexpr
real scal_prod( const point& lhs, const point& rhs ) noexcept
{
	return lhs.x*rhs.x +
	       lhs.y*rhs.y +
	       lhs.z*rhs.z;
}

/*!
 * \brief Returns the cross product of the given points.
 */
constexpr
point cross_prod( const point& lhs, const point& rhs ) noexcept
{
	return point{ lhs.y*rhs.z - lhs.z*rhs.y,
	              lhs.z*rhs.x - lhs.x*rhs.z,
	              lhs.x*rhs.y - lhs.y*rhs.x };
}

/*!
 * \brief Prints out the cartesian coordinates of the point.
 */
inline
std::ostream& operator<<( std::ostream& os, const point& rhs )
{
    return os << rhs.x << ' ' << rhs.y << ' ' << rhs.z << ' ';
}

constexpr
real triple_prod( const point&  p1, const point&  p2, const point &p3 ) noexcept
{
	return ( p1.y*p2.z - p1.z*p2.y )*p3.x +
	       ( p1.z*p2.x - p1.x*p2.z )*p3.y +
	       ( p1.x*p2.y - p1.y*p2.x )*p3.z;
}

constexpr
point grad( const tl_real &x )
{
    return point { x.derivatives[0], x.derivatives[1], x.derivatives[2] };
}

}


namespace std
{

template <>
struct hash<geometry::point>
{
    size_t operator()( const geometry::point& p ) const noexcept
    {
        using boost::hash_combine;

        std::hash<decltype(p.x)> hasher;

        size_t result { 0 };
        hash_combine( result, hasher( p.x ) );
        hash_combine( result, hasher( p.y ) );
        hash_combine( result, hasher( p.z ) );
        return result;
    }
};

}

