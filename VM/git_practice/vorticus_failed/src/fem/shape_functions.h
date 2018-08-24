/*
 * Copyright (C) 2016 Matthias Kirchhart
 *
 * This file is part of vorticus.
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
#ifndef FEM_SHAPE_FUNCTIONS_H
#define FEM_SHAPE_FUNCTIONS_H

#include "geometry/bary.h"
#include "geometry/point.h"

#include <array>

namespace fem
{

using namespace geometry;

namespace shapefcts1d
{

template <size_t order> constexpr
size_t num();

template <size_t order>
using     vals = std::array<real,num<order>()>;

template <size_t order, typename T>
using   coeffs = std::array<T,num<order>()>;

template <size_t order> constexpr
coeffs<order,real> lattice();

template <size_t order> constexpr
vals<order>  N( real xi );

template <size_t order> constexpr
vals<order> dNdxi( real xi );

template <size_t order, typename T>
T apply( const coeffs<order,T> &c, const vals<order>& v );

}

namespace shapefcts2d
{

template <size_t order>
constexpr size_t num();

template <size_t order>
using     vals = std::array<real,num<order>()>;

template <size_t order, typename T>
using   coeffs = std::array<T,num<order>()>;

template <size_t order> constexpr
coeffs<order,bary2d> lattice();

template <size_t order> constexpr
vals<order> N( bary2d p );

template <size_t order> constexpr
vals<order> dNdxi( bary2d p );

template <size_t order> constexpr
vals<order> dNdeta( bary2d p );

template <size_t order, typename T>
T apply( const coeffs<order,T> &c, const vals<order>& v );

}

namespace shapefcts3d
{

template <size_t order>
constexpr size_t num();

template <size_t order>
using     vals = std::array<real,num<order>()>;

template <size_t order, typename T>
using   coeffs = std::array<T,num<order>()>;


template <size_t order> constexpr
coeffs<order,bary3d> lattice();

template <size_t order> constexpr
vals<order> N( bary3d p );


template <size_t order> constexpr
vals<order> dNdxi( bary3d p );

template <size_t order> constexpr
vals<order> dNdeta( bary3d p );

template <size_t order> constexpr
vals<order> dNdzeta( bary3d p );

template <size_t order> constexpr
coeffs<order,point> dN( bary3d p );

template <size_t order, typename T>
T apply( const coeffs<order,T> &c, const vals<order>& v );

}

}

#include "fem/shape_functions.tpp"
#endif

