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
#ifndef BEM_SPARSE_MATRIX_H
#define BEM_SPARSE_MATRIX_H

#include "bem/multigrid.h"
#include "bem/grid_function.h"

#include <unordered_map>
#include <boost/functional/hash.hpp>

namespace bem
{


template <uint gorder, uint aorder>
class sparse_matrix
{
public:
    real& operator()( point p_row, point p_col );
    real  operator()( point p_row, point p_col ) const noexcept;
    real    get_elem( point p_row, point p_col ) const noexcept;

    template <typename T>
    auto apply( const grid_function<T,gorder,aorder> &f ) const
    -> grid_function<T,gorder,aorder>;

    template <typename T>
    auto apply_inverse_diag( const grid_function<T,gorder,aorder> &f ) const
    -> grid_function<T,gorder,aorder>;

private:
    using idx_t = std::pair<point,point>;
    struct idx_hash
    {
        size_t operator()( idx_t pp ) const noexcept
        {
            using boost::hash_combine;

            std::hash<point> hasher;
            size_t result { 0 };
            hash_combine( result, hasher(pp.first)  );
            hash_combine( result, hasher(pp.second) );
            return result;
        }
    };

    std::unordered_map< idx_t, real, idx_hash > entries_;
};

}

#include "bem/sparse_matrix.tpp"
#endif

