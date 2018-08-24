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
#ifndef BEM_DIFFERENTIALS_H
#define BEM_DIFFERENTIALS_H

#include "bem/multigrid.h"
#include "bem/shape_functions.h"

namespace bem
{

template <uint gorder, uint aorder>
auto  rotN( const triangle<gorder> &t, const bary pos ) noexcept
-> shapefcts2d::coeffs<aorder,point>;

template <uint gorder, uint aorder>
auto  rotNvec( const triangle<gorder> &t, const bary pos ) noexcept
-> shapefcts2d::coeffs<aorder,std::array<point,3>>;

}

#include "bem/differentials.tpp"
#endif

