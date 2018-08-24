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

namespace bem
{

template <uint gorder, uint aorder> inline
real& sparse_matrix<gorder,aorder>::operator()( point p_row, point p_col )
{
    return entries_[ idx_t(p_row,p_col) ];
}

template <uint gorder, uint aorder> inline
real sparse_matrix<gorder,aorder>::operator()( point p_row, point p_col ) const noexcept
{
    const auto iter = entries_.find( idx_t(p_row,p_col) );

    if ( iter == entries_.end() ) return 0;
    else return iter->second;
}

template <uint gorder, uint aorder> inline
real sparse_matrix<gorder,aorder>::get_elem( point p_row, point p_col ) const noexcept
{
    return (*this)( p_row, p_col );
}

template <uint gorder, uint aorder>
template <typename T>
auto sparse_matrix<gorder,aorder>::apply( const grid_function<T,gorder,aorder> &f ) const
-> grid_function<T,gorder,aorder>
{
    grid_function<T,gorder,aorder> result = f;
    for ( auto& value_pair: result )
        value_pair.second = T {};

    for ( const auto& entry: entries_ )
    {
        const point target_dof = entry.first.first;
        const point source_dof = entry.first.second;
        const real  value      = entry.second;

        result( target_dof ) += value*f( source_dof );
    }

    return std::move( result );
}

template <uint gorder, uint aorder>
template <typename T>
auto sparse_matrix<gorder,aorder>::apply_inverse_diag( const grid_function<T,gorder,aorder> &f ) const
-> grid_function<T,gorder,aorder>
{
    grid_function<T,gorder,aorder> result = f;
    for ( auto& value_pair: result )
    {
        point pos = value_pair.first;
        value_pair.second /= get_elem(pos,pos);
    }

    return std::move( result );
}

}

