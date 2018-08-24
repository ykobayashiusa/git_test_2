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
 * details.  *
 * You should have received a copy of the GNU General Public License along with
 * vorticus; see the file COPYING.  If not see http://www.gnu.org/licenses.
 */


namespace fem
{

template <typename T, size_t gorder, size_t aorder>
grid_function<T,gorder,aorder>::grid_function( const node_numbering<gorder,aorder> &numb ):
  numb_ ( numb ),
values_ ( numb.size() )
{}

template <typename T, size_t gorder, size_t aorder>
T grid_function<T,gorder,aorder>::eval( const tetrahedron<gorder>& t, bary3d pos ) const
{
    using shapefcts3d::N;
    using shapefcts3d::apply;

    coeffs_t  coeffs {};
    const numbs_t numbs { numb_(t) };
    for ( size_t i { 0 }; i < numbs.size(); ++i )
    {
        coeffs[i] = values_[ numbs[i] ];
    }
    return apply<aorder,T>( coeffs, N<aorder>(pos) );
}

template <typename T, size_t gorder, size_t aorder> inline
T grid_function<T,gorder,aorder>::eval( point pos ) const
{
    return values_[ numb_(pos) ];
}

template <typename T, size_t gorder, size_t aorder> inline
      T& grid_function<T,gorder,aorder>::operator()( point pos )
{
    return values_[ numb_(pos) ];
}

template <typename T, size_t gorder, size_t aorder> inline
const T& grid_function<T,gorder,aorder>::operator()( point pos ) const
{
    return values_[ numb_(pos) ];
}

template <typename T, size_t gorder, size_t aorder> inline
void grid_function<T,gorder,aorder>::operator+=( const grid_function& rhs )
{
    for ( auto& val: values_ )
        val.second += rhs( val.first );
}

template <typename T, size_t gorder, size_t aorder> inline
void grid_function<T,gorder,aorder>::operator-=( const grid_function& rhs )
{
    for ( auto& val: values_ )
        val.second -= rhs( val.first );
}

template <typename T, size_t gorder, size_t aorder> inline
void grid_function<T,gorder,aorder>::operator*=( const real rhs )
{
    for ( T& val: values_ )
        val *= rhs;
}

template <typename T, size_t gorder, size_t aorder> inline
void grid_function<T,gorder,aorder>::operator/=( const real rhs )
{
    const real q { 1/rhs };
    for ( T& val: values_ )
        val *= q;
}

template <typename T, size_t gorder, size_t aorder> inline
auto grid_function<T,gorder,aorder>::operator=( const T& rhs ) -> grid_function&
{
    for ( T& val: values_ )
        val = rhs;
    return *this;
}

}

