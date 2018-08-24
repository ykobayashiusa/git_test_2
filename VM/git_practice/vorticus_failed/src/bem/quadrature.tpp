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
#include <cmath>
#include <utility>
#include <iostream>

#include "bem/legendre.h"
#include "bem/dunavant.h"

namespace bem
{

///////////////////////////////////////////////////////////
// Declaration of helper functions and data structures.  //
///////////////////////////////////////////////////////////

constexpr real C5 {0.25};
constexpr real C6 {0.25};

struct intersection
{
    static constexpr char none  = 0;
    static constexpr char left  = 1;
    static constexpr char right = 2;

    static constexpr char disjoint = 0;
    static constexpr char corner   = 1;
    static constexpr char edge     = 2;
    static constexpr char equal    = 3;

    char type { disjoint };
    char rot1 { none };
    char rot2 { none };
};

template <uint gorder>
intersection determine_rotation( const triangle<gorder>& t1, const triangle<gorder>& t2 );




template <uint gorder, typename fct_t>
auto integrate_disjoint( fct_t F, const triangle<gorder>& t1, const triangle<gorder>& t2 )
     -> decltype( F( t1, bary {}, t2, bary {} ) );

template <uint gorder, typename fct_t>
auto integrate_singular ( fct_t F, const triangle<gorder>& t1, const triangle<gorder>& t2, intersection inter )
     -> decltype( F( t1, bary {}, t2, bary {} ) );




template <uint gorder, typename fct_t>
auto vertex_integrand( fct_t F,
                       const triangle<gorder>& t1, const triangle<gorder>& t2,
                       real eta1, real eta2, real eta3, real xi, intersection inter )
     -> decltype( F( t1, bary {}, t2, bary {} ) );
 
template <uint gorder, typename fct_t>
auto edge_integrand( fct_t F,
                     const triangle<gorder>& t1, const triangle<gorder>& t2,
                     real eta1, real eta2, real eta3, real xi, intersection inter )
     -> decltype( F( t1, bary {}, t2, bary {} ) );
 
template <uint gorder, typename fct_t>
auto self_integrand( fct_t F,
                     const triangle<gorder>& t1, const triangle<gorder>& t2,
                     real eta1, real eta2, real eta3, real xi, intersection inter )
     -> decltype( F( t1, bary {}, t2, bary {} ) );
 
template <uint gorder, typename fct_t>
auto K_3( fct_t F, const triangle<gorder>& t1, const triangle<gorder>& t2,
          real eta1, real eta2, real eta3, real xi, intersection inter )
     -> decltype( F( t1, bary {}, t2, bary {} ) );

/////////////////////
// Implementation. //
/////////////////////

template <uint gorder>
intersection determine_rotation( const triangle<gorder>& t1, const triangle<gorder>& t2 )
{
    intersection res;

    // Check for equality.
    if ( &t1 == &t2 )
    {
        res.type = intersection::equal;
        return res;
    }

    // Check for shared edges.
    for ( unsigned char e_t1 { 0 }; e_t1 < 3; ++e_t1 )
    {
        for ( unsigned char e_t2 { 0 }; e_t2 < 3; ++e_t2 )
        {
            if ( t1.get_edge_id( e_t1 ) == t2.get_edge_id( e_t2 ) )
            {
                     if ( e_t1 == 0 ) res.rot1 = intersection::left;
                else if ( e_t1 == 1 ) res.rot1 = intersection::right;
                ;
                     if ( e_t2 == 0 ) res.rot2 = intersection::left;
                else if ( e_t2 == 1 ) res.rot2 = intersection::right;
            
                res.type = intersection::edge;
                return res;
            }
        }
    }

    // Check for shared corners.
    for ( unsigned char c_t1 { 0 }; c_t1 < 3; ++c_t1 )
    {
        for ( unsigned char c_t2 { 0 }; c_t2 < 3; ++c_t2 )
        {
            if ( t1.get_node( c_t1 ) == t2.get_node( c_t2 ) )
            {
                     if ( c_t1 == 1 ) res.rot1 = intersection::left;
                else if ( c_t1 == 2 ) res.rot1 = intersection::right;
                ;
                     if ( c_t2 == 1 ) res.rot2 = intersection::left;
                else if ( c_t2 == 2 ) res.rot2 = intersection::right;

                res.type = intersection::corner;
                return res;
            }
        }
    }

    return res;
}

template <uint gorder, typename fct_t>
auto double_integral( fct_t F, const triangle<gorder>& t1, const triangle<gorder>& t2 )
     -> decltype( F( t1, bary {}, t2, bary {} ) )
{
    auto inter = determine_rotation<gorder>( t1, t2 );

    if ( inter.type == intersection::disjoint )
    {
        return integrate_disjoint( F, t1, t2 );
    }
    else
    {
        return integrate_singular( F, t1, t2, inter );
    }
}

template <uint gorder, typename fct_t>
auto integrate_disjoint   ( fct_t F, const triangle<gorder>& t1, const triangle<gorder>& t2 )
     -> decltype( F( t1, bary {}, t2, bary {} ) )
{
    using std::min;
    using std::max;
    using std::ceil;

    const real          h = max( t1.diam(), t2.diam() );
    const real          d = dist( t1, t2 );
    const unsigned char J = max( t1.level(), t2.level() );

    // The multiplication by 2 and substraction of is to convert the number of
    // Gaussian quadrature nodes to the corresponding degree of exactness.
    uint o = ( d < 3*h ) ? static_cast<uint>( ceil(C6*(J+2))*2 - 1 ) : 2u*2u - 1u;
         o = min( o, 20u );
         o = max( o, 2u*2u - 1u );

    // The following three lines create a zero-value. Necessary for math::matrix,
    // as the default constructor doesn’t know the correct dimensionse.
    auto result = F( t1, bary { 1./3., 1./3., 1./3. },
                     t2, bary { 1./3., 1./3., 1./3. } ); 
    result -= result;

    for ( uint i = 0; i < dunavant_num_points[ o ]; ++i )
    {
        real  xi1 = dunavant_points [ o ][ i ].xi();
        real eta1 = dunavant_points [ o ][ i ].eta();
        real   w1 = dunavant_weights[ o ][ i ];

        for ( uint j = 0; j < dunavant_num_points[ o ]; ++j )
        {
            real  xi2 = dunavant_points [ o ][ j ].xi();
            real eta2 = dunavant_points [ o ][ j ].eta();
            real   w2 = dunavant_weights[ o ][ j ];

            result += 0.25*w1*w2*K_3( F, t1, t2, xi1, eta1, xi2, eta2, intersection {} );
        }
    }
    return result;
}

template <uint gorder, typename fct_t>
auto integrate_singular( fct_t F,
                         const triangle<gorder>& t1, const triangle<gorder>& t2, intersection inter )
     -> decltype( F( t1, bary {}, t2, bary {} ) )
{
    using std::min;
    using std::max;
    using std::ceil;

    // The following three lines create a zero-value. Necessary for math::matrix,
    // as the default constructor doesn’t know the correct dimensionse.
    // The points are chosen differently in the inside of the triangles in order
    // to avoid hitting the singularity.
    auto result = F( t1, bary { 2./6., 3./6., 1./6. },
                     t2, bary { 2./6., 1./6., 3./6. } ); 
    result -= result;

    const unsigned char J = max( t1.level(), t2.level() );

    const uint n1 = 4u;
          uint n2 = ceil(2*C5*(J+2));
               n2 = min( n2, 20u );
               n2 = max( n2, 4u );

    for ( uint i = 0; i < n1; ++i )
    {
        real xi = legendre_points [ n1 ][ i ];
        real w1 = legendre_weights[ n1 ][ i ];
        for ( uint j = 0; j < n2; ++j )
        {
            real eta3 = legendre_points [ n2 ][ j ];
            real w2   = legendre_weights[ n2 ][ j ];
            for ( uint k = 0; k < n2; ++k )
            {
                real eta2 = legendre_points [ n2 ][ k ];
                real w3   = legendre_weights[ n2 ][ k ];
                for ( uint l = 0; l < n2; ++l )
                {
                    real eta1 = legendre_points [ n2 ][ l ];
                    real w4   = legendre_weights[ n2 ][ l ];
                   
                         if ( inter.type == intersection::corner )  result += w1 * w2 * w3 * w4 * vertex_integrand( F, t1, t2, eta1, eta2, eta3, xi, inter );
                    else if ( inter.type == intersection::edge   )  result += w1 * w2 * w3 * w4 *   edge_integrand( F, t1, t2, eta1, eta2, eta3, xi, inter );
                    else if ( inter.type == intersection::equal  )  result += w1 * w2 * w3 * w4 *   self_integrand( F, t1, t2, eta1, eta2, eta3, xi, inter );
                }
            }
        }
    }
    return result;
}

template <uint gorder, typename fct_t> inline
auto self_integrand( fct_t F,
                     const triangle<gorder>& t1, const triangle<gorder>& t2,
                     real eta1, real eta2, real eta3, real xi, intersection inter )
     -> decltype( F( t1, bary {}, t2, bary {} ) )
{
    return xi*xi*xi*eta1*eta1*eta2*
           (
              K_3( F, t1, t2, xi, xi*(1-eta1+eta1*eta2), xi*(1-eta1*eta2*eta3), xi*(1-eta1), inter )      +
              K_3( F, t1, t2, xi*(1 - eta1*eta2*eta3), xi*(1-eta1), xi, xi*(1-eta1+eta1*eta2), inter )    +
              K_3( F, t1, t2, xi, xi*eta1*(1-eta2+eta2*eta3), xi*(1-eta1*eta2), xi*eta1*(1-eta2), inter ) +
              K_3( F, t1, t2, xi*(1-eta1*eta2), xi*eta1*(1-eta2), xi, xi*eta1*(1-eta2+eta2*eta3), inter ) + 
              K_3( F, t1, t2, xi*(1-eta1*eta2*eta3), xi*eta1*(1-eta2*eta3), xi, xi*eta1*(1-eta2), inter ) +
              K_3( F, t1, t2, xi, xi*eta1*(1-eta2), xi*(1-eta1*eta2*eta3), xi*eta1*(1-eta2*eta3), inter )
           );
}

template <uint gorder, typename fct_t> inline
auto edge_integrand( fct_t F,
                     const triangle<gorder>& t1, const triangle<gorder>& t2,
                     real eta1, real eta2, real eta3, real xi, intersection inter )
     -> decltype( F( t1, bary {}, t2, bary {} ) )
{
    return xi*xi*xi*eta1*eta1*K_3( F, t1, t2, xi, xi*eta1*eta3, xi*(1-eta1*eta2), xi*eta1*(1-eta2), inter ) +
           xi*xi*xi*eta1*eta1*eta2*
           (
                K_3( F, t1, t2, xi, xi*eta1, xi*(1-eta1*eta2*eta3), xi*eta1*eta2*(1-eta3), inter ) +
                K_3( F, t1, t2, xi*(1-eta1*eta2), xi*eta1*(1-eta2), xi, xi*eta1*eta2*eta3, inter ) +
                K_3( F, t1, t2, xi*(1-eta1*eta2*eta3), xi*eta1*eta2*(1-eta3), xi, xi*eta1, inter ) +
                K_3( F, t1, t2, xi*(1-eta1*eta2*eta3), xi*eta1*(1-eta2*eta3), xi, xi*eta1*eta2, inter )
           );
}

template <uint gorder, typename fct_t> inline
auto vertex_integrand( fct_t F,
                       const triangle<gorder>& t1, const triangle<gorder>& t2,
                       real eta1, real eta2, real eta3, real xi, intersection inter )
     -> decltype( F( t1, bary {}, t2, bary {} ) )
{
    return xi*xi*xi*eta2*
           (
                K_3( F, t1, t2, xi, xi*eta1, xi*eta2, xi*eta2*eta3, inter ) +
                K_3( F, t1, t2, xi*eta2, xi*eta2*eta3, xi, xi*eta1, inter )
           );
}

template <uint gorder, typename fct_t> inline
auto K_3( fct_t F, const triangle<gorder>& t1, const triangle<gorder>& t2,
          real xi1, real eta1, real xi2, real eta2, intersection inter )
     -> decltype( F( t1, bary {}, t2, bary {} ) )
{
    using std::swap;
    using geometry::point;

    bary   x { xi1, eta1 };
    if ( inter.rot1 == intersection::left )
    {
        bary tmp = x;
        x.z0 = tmp.z2;
        x.z1 = tmp.z0;
        x.z2 = tmp.z1;
    }
    else if ( inter.rot1 == intersection::right )
    {
        bary tmp = x;
        x.z0 = tmp.z1;
        x.z1 = tmp.z2;
        x.z2 = tmp.z0;
    }

    bary   y { xi2, eta2 };
    if ( inter.rot2 == intersection::left )
    {
        bary tmp = y;
        y.z0 = tmp.z2;
        y.z1 = tmp.z0;
        y.z2 = tmp.z1;
    }
    else if ( inter.rot2 == intersection::right )
    {
        bary tmp = y;
        y.z0 = tmp.z1;
        y.z1 = tmp.z2;
        y.z2 = tmp.z0;
    }
    
    return F( t1, x, t2, y )*t1.surf_elem( x )*t2.surf_elem( y );
}

template <uint gorder, uint aorder>
math::vector rhs( const triangle<gorder>& t )
{
    using shapefcts2d::N;

    constexpr uint o { 3u*aorder };

    math::vector result { lattice::size<aorder>(), arma::fill::zeros };

    for ( uint i = 0; i < dunavant_num_points[ o ]; ++i )
    {
        bary x = dunavant_points [ o ][ i ];
        real w = dunavant_weights[ o ][ i ];

        auto Nt = N<aorder>( x );
        real  gx { t.surf_elem( x ) };
        point n  { t.normal( x ) };

        arma::colvec Nvec( &Nt[ 0 ], Nt.size(), false, true );

        point u_infty { 1, 0, 0 };
        result += 0.5*w*gx*( scal_prod( -u_infty, n ) )*Nvec;
    }
    return result;
}

template <uint gorder, uint aorder>
auto elem_mass( const triangle<gorder>& t )
-> arma::mat::fixed< shapefcts2d::num<aorder>(), shapefcts2d::num<aorder>() >
{
    using shapefcts2d::N;
    using shapefcts2d::num;
    constexpr size_t size = num<aorder>();

    arma::mat::fixed<size,size> result;
    result.zeros();

    constexpr uint o = (aorder == 0u) ? 3u : 3u*aorder;
    
    for ( uint i = 0; i < dunavant_num_points[ o ]; ++i )
    {
        const bary pos = dunavant_points [ o ][ i ];
        const real w   = dunavant_weights[ o ][ i ];
        const real gx  = t.surf_elem( pos );

        auto NN = N<aorder>( pos );
        const arma::Col<real> NN1( &(NN[0]), size, false, true );
        const arma::Row<real> NN2( &(NN[0]), size, false, true );
        arma::mat::fixed<size,size> tmp = NN1*NN2;
        tmp *= 0.5*w*gx;

        result += tmp;
    }

    return result;
}

template <uint gorder, uint aorder>
sparse_matrix<gorder,aorder> mass_matrix( const grid<gorder>& g )
{
    grid_function<real,gorder,aorder> f(g);
    sparse_matrix<gorder,aorder> result;

    for ( const triangle<gorder>& t: g )
    {
        const auto& dofs = f.get_dofs( t );
        const auto  elem_mat = elem_mass<gorder,aorder>( t );

        for ( size_t i = 0; i < dofs.size(); ++i )
        {
            for ( size_t j = 0; j < dofs.size(); ++j )
            {
                result( dofs[i], dofs[j] ) += elem_mat(i,j);
            }
        }
    }
    return std::move(result);
}

}

