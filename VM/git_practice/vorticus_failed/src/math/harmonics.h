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
#ifndef MATH_HARMONICS_H
#define MATH_HARMONICS_H

#include "types.h"

#include "geometry/point.h"
#include "geometry/cmplx_point.h"

#include <cmath>
#include <array>

namespace math
{

///////////////////////
// Type definitions. //
///////////////////////

template <int order>
struct real_coeff
{
    static constexpr int num_coeff();
        
          real* begin()       { return data.begin(); }
          real* end()         { return data.end();   }
    const real* begin() const { return data.begin(); }
    const real* end()   const { return data.end();   }

    real operator()( int i, int j ) const;

    std::array<real,num_coeff()> data {{}};
};

template <int order>
struct cmplx_coeff
{
    static constexpr int num_coeff();
        
          cmplx* begin()       { return data.begin(); }
          cmplx* end()         { return data.end();   }
    const cmplx* begin() const { return data.begin(); }
    const cmplx* end()   const { return data.end();   }

    cmplx operator()( int i, int j ) const;

    void operator+=( const cmplx_coeff<order> &rhs );

    std::array<cmplx,num_coeff()> data {{}};
};



template <int order>
struct cmplx_point_coeff
{
    static constexpr int num_coeff();

          geometry::cmplx_point* begin()       { return data.begin(); }
          geometry::cmplx_point* end()         { return data.end();   }
    const geometry::cmplx_point* begin() const { return data.begin(); }
    const geometry::cmplx_point* end()   const { return data.end();   }

    geometry::cmplx_point operator()( int i, int j ) const;

    std::array<geometry::cmplx_point,num_coeff()> data {{}};
};



//////////////////////////////////
// Declarations of functions.   //
//////////////////////////////////

constexpr
int idx( const int n, const int m );

template <int order>
cmplx_coeff<order> R( const geometry::point r );

template<int order>
cmplx_point_coeff<order> dR( const geometry::point r );

template <int order>
cmplx_coeff<order> S( const geometry::point r );

template<int order>
cmplx_point_coeff<order> dS( const geometry::point r );

/////////////////////////////////////
// Member function implementation. //
/////////////////////////////////////

template <int order>
constexpr int real_coeff<order>::num_coeff()
{
    return (order*(order+1))/2;
}

template <int order> inline
real real_coeff<order>::operator()( int n, int m ) const
{
    return data[ idx(n,m) ];
}

template <int order> inline
void cmplx_coeff<order>::operator+=( const cmplx_coeff<order> &rhs )
{
    for ( int i = 0; i < num_coeff(); ++i )
        data[ i ] += rhs.data[ i ];
}

template <int order>
constexpr int cmplx_coeff<order>::num_coeff()
{
    return (order*(order+1))/2;
}

template <int order> inline
cmplx cmplx_coeff<order>::operator()( int n, int m ) const
{
    if ( m < 0 )
    {
        if ( m % 2 ) return -conj( data[ idx(n,-m) ] );
        else         return  conj( data[ idx(n,-m) ] );
    }
    else return data[ idx(n,m) ];
}

template <int order>
constexpr int cmplx_point_coeff<order>::num_coeff()
{
    return (order*(order+1))/2;
}

template <int order> inline
geometry::cmplx_point cmplx_point_coeff<order>::operator()( int n, int m ) const
{
    if ( m < 0 )
    {
        if ( m % 2 ) return -conj( data[ idx(n,-m) ] );
        else         return  conj( data[ idx(n,-m) ] );
    }
    else return data[ idx(n,m) ];
}

/////////////////////////////////////////
// Non-member function implementation. //
/////////////////////////////////////////

constexpr
int idx( const int n, const int m )
{
    return m + ( n*(n+1) )/2;
}

template <int order>
cmplx_coeff<order> R( const geometry::point p_r )
{
    const real  x = p_r.x;
    const real  y = p_r.y;
    const real  z = p_r.z;
    const real r2 = x*x + y*y + z*z;
    const cmplx xy { x, y };

    cmplx_coeff<order> result;
    result.data[ idx(0,0) ] = cmplx { 1,0 };

    // First do the diagonal.
    for ( int m = 0; m < order - 1; ++m )
    {
        result.data[ idx(m+1,m+1) ] = -result.data[ idx(m,m) ]*xy/((real)(2*m+2));
    }

    // Afterwards the "small diagonal" below.
    for ( int m = 0; m < order - 1; ++m )
    {
        result.data[ idx(m+1,m) ] = result.data[ idx(m,m) ]*z;
    }

    // Now the rest.
    for ( int m = 0; m < order - 2; ++m )
    {
        for ( int n = m + 2; n < order; ++n )
        {
            result.data[ idx(n,m) ] =
            (
                result.data[ idx(n-1,m) ]*(z*(2*n-1)) -
                result.data[ idx(n-2,m) ]*r2
            )/((real)((n+m)*(n-m)));
        }
    }

    return result;
}

template<int order>
cmplx_point_coeff<order> dR( const geometry::point p_r )
{
    cmplx_coeff<order> r = R<order>( p_r );

    cmplx_point_coeff<order> result;

    // First do the diagonal.
    for ( int m = 0; m < order - 1; ++m )
    {
        result.data[ idx(m+1,m+1) ].x = -r.data[ idx(m,m) ]/((real) 2);
        result.data[ idx(m+1,m+1) ].y.real( r.data[ idx(m,m) ].imag() /((real) 2) );
        result.data[ idx(m+1,m+1) ].y.imag(-r.data[ idx(m,m) ].real() /((real) 2) );
        result.data[ idx(m+1,m+1) ].z = 0;
    }

    // Afterwards the "small diagonal" below.
    result.data[ idx(1,0) ].x = 0;
    result.data[ idx(1,0) ].y = 0;
    result.data[ idx(1,0) ].z = 1;
    for ( int m = 1; m < order - 1; ++m )
    {
        result.data[ idx(m+1,m) ].x = -r.data[ idx(m,m-1) ]/((real) 2);
        result.data[ idx(m+1,m) ].y.real( r.data[ idx(m,m-1) ].imag() /((real) 2) );
        result.data[ idx(m+1,m) ].y.imag(-r.data[ idx(m,m-1) ].real() /((real) 2) );
        result.data[ idx(m+1,m) ].z = r.data[ idx(m,m) ];
    }

    // Now the rest.
    for ( int m = 0; m < order - 2; ++m )
    {
        for ( int n = m + 2; n < order; ++n )
        {
            result.data[ idx(n,m) ].x = ( r(n-1,m+1) - r(n-1,m-1) ) / ((real) 2);
            result.data[ idx(n,m) ].y.real( (r(n-1,m+1).imag() + r(n-1,m-1).imag()) / ((real) 2) );
            result.data[ idx(n,m) ].y.imag(-(r(n-1,m+1).real() + r(n-1,m-1).real()) / ((real) 2) );
            result.data[ idx(n,m) ].z = r(n-1,m);
        }
    }

    return result;
}

template <int order>
cmplx_coeff<order> S( const geometry::point p_r )
{
    const real  x = p_r.x;
    const real  y = p_r.y;
    const real  z = p_r.z;
    const real r2 = x*x + y*y + z*z;
    const cmplx xy { x, y };

    cmplx_coeff<order> result;
    result.data[ idx(0,0) ] = cmplx { 1/std::sqrt(r2),0 };

    // First do the diagonal.
    for ( int m = 0; m < order - 1; ++m )
    {
        result.data[ idx(m+1,m+1) ] = -((2*m+1)/r2)*xy*result.data[ idx(m,m) ];
    }

    // Afterwards the "small diagonal" below.
    for ( int m = 0; m < order - 1; ++m )
    {
        result.data[ idx(m+1,m) ] = ((2*m+1)/r2)*result.data[ idx(m,m) ]*z;
    }

    // Now the rest.
    for ( int m = 0; m < order - 2; ++m )
    {
        for ( int n = m + 2; n < order; ++n )
        {
            result.data[ idx(n,m) ] = 
            (
                result.data[ idx(n-1,m) ]*(z/r2)*((real)(2*n-1)) -
                result.data[ idx(n-2,m) ]*(((n+m-1)*(n-m-1))/r2)
            );
        }
    }

    return result;
}

template<int order>
cmplx_point_coeff<order> dS( const geometry::point p_r )
{
    cmplx_coeff<order+1> s = S<order+1>( p_r );

    cmplx_point_coeff<order> result;
    for ( int n = 0; n < order; ++n )
    {
        for ( int m = 0; m <= n; ++m )
        {
            result.data[ idx(n,m) ].x = ( s(n+1,m+1) - s(n+1,m-1) )/((real) 2);
            result.data[ idx(n,m) ].y.real( (s(n+1,m+1).imag() + s(n+1,m-1).imag())/((real) 2) );
            result.data[ idx(n,m) ].y.imag(-(s(n+1,m+1).real() + s(n+1,m-1).real())/((real) 2) );
            result.data[ idx(n,m) ].z = -s(n+1,m);
        }
    }

    return result;
}

} // namespace math.

#endif

