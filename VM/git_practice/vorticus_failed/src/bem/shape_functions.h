/*
 * Copyright (C) 2014 Matthias Kirchhart
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
#ifndef BEM_SHAPE_FUNCTIONS_H
#define BEM_SHAPE_FUNCTIONS_H

#include "types.h"
#include "bary.h"

#include <array>

namespace bem
{

namespace shapefcts1d
{

template <uint order> constexpr
uint num();

template <uint order>
using     vals = std::array<real,num<order>()>;

template <uint order, typename T>
using   coeffs = std::array<T,num<order>()>;

template <uint order> constexpr
vals<order>  N( real xi );

template <uint order> constexpr
vals<order> dN( real xi );

template <uint order, typename T>
T apply( const coeffs<order,T> &c, const vals<order>& v );

}

namespace shapefcts2d
{

template <uint order>
constexpr uint num();

template <uint order>
using     vals = std::array<real,num<order>()>;

template <uint order, typename T>
using   coeffs = std::array<T,num<order>()>;

template <uint order> constexpr
vals<order> N( bary p );

template <uint order> constexpr
vals<order> dNdxi( bary p );

template <uint order> constexpr
vals<order> dNdeta( bary p );

template <uint order, typename T>
T apply( const coeffs<order,T> &c, const vals<order>& v );

}

}

#include "shape_functions.tpp"
#endif

