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
 * \file geometry/point.h
 * \author Matthias Kirchhart <matthias@kirchhart.de>
 * \brief Definition of the geometry::point class.
 */
#ifndef GEOMETRY_POINT_H
#define GEOMETRY_POINT_H

#include <cmath>
#include <cfloat>
#include <ostream>
#include <functional>
#include <boost/functional/hash.hpp>

#include "types.h"

namespace geometry
{

/*!
 * \brief A three-dimensional vector in Euclidian space.
 *
 * Objects of this class represent a vector in the three-dimensional Euclidean
 * space. They can also thought of vectors \f$\vec{v}\in\mathbb{R}^3\f$.
 * The coordinates are stored in Cartesian format, but methods for conversion to
 * spherical coordinates also exist. This class is called <i>point</i> in order
 * to emphasize its geometric nature. For general tuples of numbers, the
 * math::vector class should be used.
 */
struct point 
{
	constexpr point() noexcept = default;
	constexpr point( real xx, real yy, real zz ) noexcept;

	real       r () const noexcept;
	real    theta() const noexcept;
	real      phi() const noexcept;

	constexpr real  r2() const noexcept;

    real max_norm() const noexcept;

	point e_r()     const noexcept;
	point e_theta() const noexcept;
	point e_phi()   const noexcept;

	constexpr point operator-() const noexcept;
	
	constexpr point operator+( const point& rhs ) const noexcept;
	constexpr point operator-( const point& rhs ) const noexcept;
	
	void operator+=( const point& rhs ) noexcept;
	void operator-=( const point& rhs ) noexcept;

	constexpr point operator*( real rhs ) const noexcept;
	constexpr point operator/( real rhs ) const noexcept;

	void operator*=( real rhs ) noexcept;
	void operator/=( real rhs ) noexcept;

	real x { 0 };
	real y { 0 };
	real z { 0 };
};

real length( const point& p ) noexcept;

point min_components( const point& lhs, const point& rhs ) noexcept;
point max_components( const point& lhs, const point& rhs ) noexcept;

constexpr bool operator==( const point& lhs, const point& rhs ) noexcept;
constexpr bool operator!=( const point& lhs, const point& rhs ) noexcept;

constexpr point   operator*(         real lhs, const point& rhs ) noexcept;
constexpr real    scal_prod( const point& lhs, const point& rhs ) noexcept;
constexpr point  cross_prod( const point& lhs, const point& rhs ) noexcept;
constexpr real  triple_prod( const point&  p1, const point&  p2, const point &p3 ) noexcept;

std::ostream& operator<<( std::ostream& os, const point& rhs );

constexpr point grad( const tl_real &x );

}

#include "geometry/point.tpp"
#endif

