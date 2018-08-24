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
#ifndef VRM_LAPLACIAN_H
#define VRM_LAPLACIAN_H

#include <mutex>
#include <vector>
#include <utility>
#include <algorithm>
#include <armadillo>

#include "math/simplex.h"
#include "math/rev_simplex.h"

namespace vrm
{

template <typename iterator>
using matrix_entry = std::tuple<iterator,iterator,real>;

template <typename iterator>
using stiffness_matrix = std::vector<matrix_entry<iterator>>;

template <typename iterator, typename accessor,
          typename finder,   typename maker>
auto laplacian_2d( iterator begin, iterator end,
                   accessor &get, finder &find, maker &make,
                   real h, real C_diff )
-> stiffness_matrix<iterator>;

}

#include "vrm/laplacian.tpp"
#endif

