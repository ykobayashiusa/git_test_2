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
#ifndef FEM_GRID_FUNCTION_H
#define FEM_GRID_FUNCTION_H

#include "fem/multigrid.h"
#include "fem/node_numbering.h"

namespace fem
{

template <typename T, size_t gorder>
class grid_function_base
{
public:
    virtual T eval( point pos ) const = 0;
    virtual T eval( const tetrahedron<gorder>& t, bary3d pos ) const = 0;
};


template <typename T, size_t gorder, size_t aorder>
class grid_function: public grid_function_base<T,gorder>
{
public:
    using       coords_t = shapefcts3d::coeffs<aorder,point>;
    using        numbs_t = shapefcts3d::coeffs<aorder,size_t>;
    using       coeffs_t = shapefcts3d::coeffs<aorder,T>;

    grid_function() = delete;
    grid_function( const node_numbering<gorder,aorder> &numb );
    grid_function( const grid_function&  rhs ) = default;
    grid_function(       grid_function&& rhs ) = default;
    grid_function& operator=( const grid_function&  rhs ) = default;
    grid_function& operator=(       grid_function&& rhs ) = default;
    grid_function& operator=( const T& rhs );

          T& operator()( point pos );
    const T& operator()( point pos ) const;

    T eval( point pos ) const override;
    T eval( const tetrahedron<gorder>& t, bary3d pos ) const override;


    void operator+=( const grid_function& rhs );
    void operator-=( const grid_function& rhs );
    void operator*=( const real rhs );
    void operator/=( const real rhs );

private:
    const node_numbering<gorder,aorder>& numb_;
    std::vector<T> values_;
};

}

#include "fem/grid_function.tpp"
#endif

