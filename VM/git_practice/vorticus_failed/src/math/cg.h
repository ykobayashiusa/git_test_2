/*
 * Copyright (C) 2016 Matthias Kirchhart
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
#ifndef MATH_CG_H
#define MATH_CG_H

#include "types.h"

namespace math
{

template <typename mat, typename vec>
real cg( mat& A, vec& b, vec &x,
         real target_residual = 1e-6, size_t max_iter = 1000, bool relative = true );

template <typename precond, typename mat, typename vec>
real cg( precond &M, mat& A, vec& b, vec &x,
         real target_residual = 1e-6, size_t max_iter = 1000, bool relative = true );
}

#include "math/cg.tpp"
#endif

