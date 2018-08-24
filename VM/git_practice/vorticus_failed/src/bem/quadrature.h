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
#ifndef BEM_QUADRATURE_H
#define BEM_QUADRATURE_H

#include "math/matrix.h"
#include "math/vector.h"

#include "bem/multigrid.h"
#include "bem/sparse_matrix.h"

namespace bem
{

template <uint gorder, typename fct_t>
auto double_integral( fct_t F, const triangle<gorder>& t1, const triangle<gorder>& t2 )
     -> decltype( F( t1, bary {}, t2, bary {} ) );

template <uint gorder, uint aorder>
auto    elem_mass( const triangle<gorder>& t )
-> arma::mat::fixed< shapefcts2d::num<aorder>(), shapefcts2d::num<aorder>() >;

template <uint gorder, uint aorder>
sparse_matrix<gorder,aorder> mass_matrix( const grid<gorder>& g );

template <uint gorder, uint aorder>
math::vector rhs( const triangle<gorder>& t );

}

#include "bem/quadrature.tpp"
#endif

