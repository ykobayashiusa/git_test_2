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
#ifndef MATH_POTENTIALS_H
#define MATH_POTENTIALS_H

#include "geometry/tl_point.h"
#include "fem/shape_functions.h"

namespace fem
{

template <size_t order>
std::pair< shapefcts2d::coeffs<order,real>, shapefcts2d::coeffs<order,real> >
layer_potentials( const std::array<point,3> X, point P );

template <>
std::pair< shapefcts2d::coeffs<0,real>, shapefcts2d::coeffs<0,real> >
layer_potentials<0>( const std::array<point,3> X, point P );

template <>
std::pair< shapefcts2d::coeffs<1,real>, shapefcts2d::coeffs<1,real> >
layer_potentials<1>( const std::array<point,3> X, point P );




template <size_t order>
std::pair< shapefcts2d::coeffs<order,tl_real>, shapefcts2d::coeffs<order,tl_real> >
layer_potentials( const std::array<point,3> X, geometry::tl_point P );

template <>
std::pair< shapefcts2d::coeffs<0,tl_real>, shapefcts2d::coeffs<0,tl_real> >
layer_potentials<0>( const std::array<point,3> X, geometry::tl_point P );

template <>
std::pair< shapefcts2d::coeffs<1,tl_real>, shapefcts2d::coeffs<1,tl_real> >
layer_potentials<1>( const std::array<point,3> X, geometry::tl_point P );





template <size_t order>
point biot_savart   ( const std::array<point,4> X, const shapefcts3d::coeffs<order,point> O, point P );

template <>
point biot_savart<0>( const std::array<point,4> X, const std::array<point,1> O, const point P );

template <>
point biot_savart<1>( const std::array<point,4> X, const std::array<point,4> O, const point P );

}

#endif

