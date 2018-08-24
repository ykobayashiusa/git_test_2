/*
 * Copyright (C) 2015 Matthias Kirchhart
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
#ifndef MATH_NNLS_H
#define MATH_NNLS_H

#include <armadillo>

namespace math
{

template <typename flt>
arma::Col<flt> nnls( const arma::Mat<flt> &A, const arma::Col<flt> &b,
                     flt tol = 1e-9 );

}

#include "math/nnls.tpp"
#endif

