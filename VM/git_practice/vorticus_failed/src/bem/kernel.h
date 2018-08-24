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
#ifndef BEM_KERNEL_H
#define BEM_KERNEL_H

#include "math/matrix.h"
#include "math/constants.h"

#include "bem/multigrid.h"
#include "bem/shape_functions.h"

namespace bem
{

template <uint gorder, uint aorder>
auto W_kernel( const triangle<gorder>& t1, const bary pos1,
               const triangle<gorder>& t2, const bary pos2 )
-> arma::mat::fixed< shapefcts2d::num<aorder>(), shapefcts2d::num<aorder>() >;

}

#include "bem/kernel.tpp"
#endif

