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
 * \file geometry/cmplx_point.h
 * \author Matthias Kirchhart <matthias@kirchhart.de>
 * \brief Definition of the cmplx_point class.
 */
#ifndef GEOMETRY_CMPLX_POINT_H
#define GEOMETRY_CMPLX_POINT_H

#include "types.h"
#include "geometry/point.h"

#include <cmath>

namespace geometry
{

/*!
 * \brief Complex vectors from \f$\mathbb{C}^3\f$.
 * \see geometry::point
 *
 * This class is a complex version of geometry::point. For more information
 * see there. It is mainly used to store the complex coefficients in the
 * fast multipole method.
 */
struct cmplx_point
{
	constexpr cmplx_point() noexcept = default;
	constexpr cmplx_point( const point& rhs ) noexcept;
	constexpr cmplx_point( const point& real, const point& imag ) noexcept;
	constexpr cmplx_point( cmplx xx, cmplx yy, cmplx zz ) noexcept;

	      point real() const noexcept;
	      point imag() const noexcept;
	cmplx_point conj() const noexcept;

	cmplx_point operator-() const noexcept;

	cmplx_point operator+( const cmplx_point& rhs ) const noexcept;
	cmplx_point operator-( const cmplx_point& rhs ) const noexcept;

	void operator+=( const cmplx_point& rhs ) noexcept;
	void operator-=( const cmplx_point& rhs ) noexcept;

	cmplx_point operator*( cmplx rhs ) const noexcept;
	cmplx_point operator/( cmplx rhs ) const noexcept;

	void operator*=( cmplx rhs ) noexcept;
	void operator/=( cmplx rhs ) noexcept;

	cmplx x { 0, 0 };
	cmplx y { 0, 0 };
	cmplx z { 0, 0 };
};

real length( const cmplx_point& p ) noexcept;

cmplx_point operator*( cmplx lhs, const cmplx_point& rhs ) noexcept;

cmplx_point conj( const cmplx_point& p ) noexcept;

cmplx        scal_prod( const cmplx_point& lhs, const cmplx_point& rhs ) noexcept;
cmplx_point cross_prod( const cmplx_point& lhs, const cmplx_point& rhs ) noexcept;

}


#include "geometry/cmplx_point.tpp"
#endif

