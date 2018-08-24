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
#ifndef GEOMETRY_TL_POINT_H
#define GEOMETRY_TL_POINT_H

#include <cmath>
#include <cfloat>

#include "types.h"
#include "geometry/point.h"
#include "geometry/tensor.h"

namespace geometry
{

struct tl_point 
{
	tl_point() noexcept = default;
	tl_point( point x ) noexcept;
	tl_point( tl_real xx, tl_real yy, tl_real zz ) noexcept;
	tl_point( point x, tensor dx ) noexcept;

	tl_real       r () const noexcept;
	tl_real       r2() const noexcept;

	tl_point operator-() const noexcept;
	
	tl_point operator+( const point& rhs ) const noexcept;
	tl_point operator-( const point& rhs ) const noexcept;
	tl_point operator+( const tl_point& rhs ) const noexcept;
	tl_point operator-( const tl_point& rhs ) const noexcept;
	
	void operator+=( const point& rhs ) noexcept;
	void operator-=( const point& rhs ) noexcept;
	void operator+=( const tl_point& rhs ) noexcept;
	void operator-=( const tl_point& rhs ) noexcept;

	tl_point operator*( real rhs ) const noexcept;
	tl_point operator/( real rhs ) const noexcept;
	tl_point operator*( tl_real rhs ) const noexcept;
	tl_point operator/( tl_real rhs ) const noexcept;
	void operator*=( real rhs ) noexcept;
	void operator/=( real rhs ) noexcept;
	void operator*=( tl_real rhs ) noexcept;
	void operator/=( tl_real rhs ) noexcept;

    point  val () const noexcept;
    tensor grad() const noexcept;
    point  curl() const noexcept;

	tl_real x {}; 
	tl_real y {};
	tl_real z {};
};

tl_real length( const tl_point& p ) noexcept;

tl_point operator+( const point& lhs, const tl_point &rhs ) noexcept;
tl_point operator-( const point& lhs, const tl_point &rhs ) noexcept;

tl_point operator*( const tl_real& lhs, const tl_point& rhs ) noexcept;
tl_point operator*( const    real  lhs, const tl_point& rhs ) noexcept;
tl_point operator*( const tl_real& lhs, const     point rhs ) noexcept;

tl_real scal_prod( const tl_point& lhs, const tl_point& rhs ) noexcept;
tl_real scal_prod( const tl_point& lhs, const    point& rhs ) noexcept;
tl_real scal_prod( const    point& lhs, const tl_point& rhs ) noexcept;

tl_point cross_prod( const tl_point& lhs, const tl_point& rhs ) noexcept;
tl_point cross_prod( const    point& lhs, const tl_point& rhs ) noexcept;
tl_point cross_prod( const tl_point& lhs, const    point& rhs ) noexcept;

}

#include "geometry/tl_point.tpp"
#endif

