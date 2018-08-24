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
#include "tensor.h"
#include <armadillo>

namespace geometry
{

real tensor::norm() const
{
	arma::Mat<real> tmp( 3, 3 );
	tmp(0,0) = data[0][0];
	tmp(0,1) = data[0][1];
	tmp(0,2) = data[0][2];

	tmp(1,0) = data[1][0];
	tmp(1,1) = data[1][1];
	tmp(1,2) = data[1][2];

	tmp(2,0) = data[2][0];
	tmp(2,1) = data[2][1];
	tmp(2,2) = data[2][2];
	return arma::norm(tmp,2);
}

}

