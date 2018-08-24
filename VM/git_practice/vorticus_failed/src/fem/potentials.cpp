/*
 * Copyright (C) 2017 Matthias Kirchhart
 *
 * This file is part of vorticus.
 *
 * vorticus is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3, or (at your option) any later
 * version.
 * * vorticus is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * vorticus; see the file COPYING.  If not see http://www.gnu.org/licenses.
 */
#include <limits>

#include "fem/potentials.h"
#include "math/constants.h"
#include "geometry/tl_point.h"

#include <boost/simd/pack.hpp>
//#include <boost/proto/transform/detail/pack.hpp>
#include <boost/simd/function/log.hpp>
//#include <boost/thread/detail/log.hpp>
#include <boost/simd/function/sqrt.hpp>
#include <boost/simd/function/load.hpp>
#include <boost/simd/function/store.hpp>
#include <boost/simd/function/atan2.hpp>

namespace fem
{

template <>
std::pair< shapefcts2d::coeffs<0,real>, shapefcts2d::coeffs<0,real> >
layer_potentials<0>( const std::array<point,3> X, point P )
{
    using std::log;
    using std::atan2;
    using std::array;
    using math::pifac;
    constexpr real eps { std::numeric_limits<real>::epsilon() }; 

    point E[3] = { X[1]-X[0], X[2]-X[1], X[0]-X[2] };
    point R[3] = { X[0]-P, X[1]-P, X[2]-P };
    point N    = cross_prod(E[0],E[1]);

    const real  e[3] = { E[0].r(), E[1].r(), E[2].r() };
    const real  r[3] = { R[0].r(), R[1].r(), R[2].r() };
    const real  n    =   N.r();
    E[0] *= 1/e[0]; E[1] *= 1/e[1]; E[2] *= 1/e[2]; N *= 1/n;

    const point RE[3] = { cross_prod(R[0],E[0]), cross_prod(R[1],E[1]), cross_prod(R[2],E[2]) };
    const real T = triple_prod(R[0],R[1],R[2]);
    const real d = scal_prod(R[0],N); // Signed distance of P to triangle's plane.

    const real Omega = 2*atan2
    (
        T, 
        r[0]*r[1]*r[2] +
        r[0]*scal_prod(R[1],R[2])+
        r[1]*scal_prod(R[0],R[2])+
        r[2]*scal_prod(R[0],R[1])
    );
  
    const real L[3] =
    { 
        log( ( ( (r[1]+e[0])*(r[1]+e[0]) - r[0]*r[0]) + eps ) /
             ( (-(r[0]-e[0])*(r[0]-e[0]) + r[1]*r[1]) + eps ) ),

        log( ( ( (r[2]+e[1])*(r[2]+e[1]) - r[1]*r[1]) + eps ) /
             ( (-(r[1]-e[1])*(r[1]-e[1]) + r[2]*r[2]) + eps ) ),

        log( ( ( (r[0]+e[2])*(r[0]+e[2]) - r[2]*r[2]) + eps ) /
             ( (-(r[2]-e[2])*(r[2]-e[2]) + r[0]*r[0]) + eps ) )
    };

    const real I  = Omega;
    const real II = scal_prod(N,RE[0]*L[0] + RE[1]*L[1] + RE[2]*L[2]) - d*Omega;
                        
    return std::pair< shapefcts2d::coeffs<0,real>, shapefcts2d::coeffs<0,real> >
    {
       {{ II*pifac }},
       {{ -I*pifac }}
    };
}

template <>
std::pair< shapefcts2d::coeffs<1,real>, shapefcts2d::coeffs<1,real> >
layer_potentials<1>( const std::array<point,3> X, point P )
{
    using std::log;
    using std::atan2;
    using std::array;
    using math::pifac;
    constexpr real eps { std::numeric_limits<real>::epsilon() }; 

    // Triangle-only data.
    point E[3] = { X[1]-X[0], X[2]-X[1], X[0]-X[2] };
    point R[3] = { X[0]-P,   X[1]-P,   X[2]-P   };
    point N    =  cross_prod(E[0],E[1]);
    const real  e[3] = { E[0].r(),  E[1].r(),  E[2].r()  };
    const real  r[3] = { R[0].r(), R[1].r(), R[2].r() };
    const real  n    =   N.r();
    E[0] *= 1/e[0]; E[1] *= 1/e[1]; E[2] *= 1/e[2]; N *= 1/n;

    const point RE[3] = { cross_prod(R[0],E[0]), cross_prod(R[1],E[1]), cross_prod(R[2],E[2]) };
    const real T = triple_prod(R[0],R[1],R[2]);
    const real d = scal_prod(R[0],N); // Signed distance of P to triangle's plane.

    const real Omega = 2*atan2
    (
        T,
        r[0]*r[1]*r[2] +
        r[0]*scal_prod(R[1],R[2])+
        r[1]*scal_prod(R[0],R[2])+
        r[2]*scal_prod(R[0],R[1])
    );
  
    const real L[3] =
    { 
         log( ( ( (r[1]+e[0])*(r[1]+e[0]) - r[0]*r[0]) + eps ) /
              ( (-(r[0]-e[0])*(r[0]-e[0]) + r[1]*r[1]) + eps ) ),

         log( ( ( (r[2]+e[1])*(r[2]+e[1]) - r[1]*r[1]) + eps ) /
              ( (-(r[1]-e[1])*(r[1]-e[1]) + r[2]*r[2]) + eps ) ),

         log( ( ( (r[0]+e[2])*(r[0]+e[2]) - r[2]*r[2]) + eps ) /
              ( (-(r[2]-e[2])*(r[2]-e[2]) + r[0]*r[0]) + eps ) )
    };

    const real M[3]
    {
       (r[0]+r[1])*(e[0]+(r[1]-r[0])*(r[1]-r[0])/e[0])/4 + RE[0].r2()*L[0]/2,
       (r[1]+r[2])*(e[1]+(r[2]-r[1])*(r[2]-r[1])/e[1])/4 + RE[1].r2()*L[1]/2,
       (r[2]+r[0])*(e[2]+(r[0]-r[2])*(r[0]-r[2])/e[2])/4 + RE[2].r2()*L[2]/2
    };

    const real  II  = scal_prod(N,RE[0]*L[0] + RE[1]*L[1] + RE[2]*L[2]) - d*Omega;
    const point III = Omega*N + cross_prod(N,E[0]*L[0]+E[1]*L[1]+E[2]*L[2]);
    const point IV  = II*N;
    const point V   = E[0]*M[0] + E[1]*M[1] + E[2]*M[2];

    point singlet = point( scal_prod(RE[1], IV), scal_prod(RE[2], IV), scal_prod(RE[0] ,IV) ) -
                    point( scal_prod( E[1],  V), scal_prod( E[2],  V), scal_prod( E[0],  V) );
    point doublet = point( scal_prod(RE[1],III), scal_prod(RE[2],III), scal_prod(RE[0],III) );
    singlet.x *= e[1]; singlet.y *= e[2]; singlet.z *= e[0];
    doublet.x *= e[1]; doublet.y *= e[2]; doublet.z *= e[0];
    singlet   *=  pifac/n;
    doublet   *= -pifac/n; 

    return std::pair< shapefcts2d::coeffs<1,real>, shapefcts2d::coeffs<1,real> >
    {
       {{ singlet.x, singlet.y, singlet.z }},
       {{ doublet.x, doublet.y, doublet.z }}
    };
}


template <>
std::pair< shapefcts2d::coeffs<0,tl_real>, shapefcts2d::coeffs<0,tl_real> >
layer_potentials<0>( const std::array<point,3> X, geometry::tl_point P )
{
    using math::pifac;
    using geometry::tl_point;

    constexpr real eps { std::numeric_limits<real>::epsilon() }; 

       point E[3] = { X[1]-X[0], X[2]-X[1], X[0]-X[2] };
    tl_point R[3] = { X[0]-P,    X[1]-P,    X[2]-P    };
       point N    = cross_prod(E[0],E[1]);

    const    real  e[3] = { E[0].r(), E[1].r(), E[2].r() };
    const tl_real  r[3] = { R[0].r(), R[1].r(), R[2].r() };
    const    real  n    =   N.r();
    E[0] *= 1/e[0]; E[1] *= 1/e[1]; E[2] *= 1/e[2]; N *= 1/n;

    const tl_point RE[3] = { cross_prod(R[0],E[0]), cross_prod(R[1],E[1]), cross_prod(R[2],E[2]) };
    const tl_real  T     = scal_prod(R[0],cross_prod(R[1],R[2]));
    const tl_real  d     = scal_prod(R[0],N); // Signed distance of P to triangle's plane.

    const tl_real Omega = 2*atan2
    (
        T, 
        r[0]*r[1]*r[2] +
        r[0]*scal_prod(R[1],R[2])+
        r[1]*scal_prod(R[0],R[2])+
        r[2]*scal_prod(R[0],R[1])
    );
  
    const tl_real L[3] =
    { 
        log( ( ( (r[1]+e[0])*(r[1]+e[0]) - r[0]*r[0]) + eps ) /
             ( (-(r[0]-e[0])*(r[0]-e[0]) + r[1]*r[1]) + eps ) ),

        log( ( ( (r[2]+e[1])*(r[2]+e[1]) - r[1]*r[1]) + eps ) /
             ( (-(r[1]-e[1])*(r[1]-e[1]) + r[2]*r[2]) + eps ) ),

        log( ( ( (r[0]+e[2])*(r[0]+e[2]) - r[2]*r[2]) + eps ) /
             ( (-(r[2]-e[2])*(r[2]-e[2]) + r[0]*r[0]) + eps ) )
    };

    const tl_real I  = Omega;
    const tl_real II = scal_prod(N,RE[0]*L[0] + RE[1]*L[1] + RE[2]*L[2]) - d*Omega;
                        
    return std::pair< shapefcts2d::coeffs<0,tl_real>, shapefcts2d::coeffs<0,tl_real> >
    {
       {{ II*pifac }},
       {{ -I*pifac }}
    };
}

template <>
std::pair< shapefcts2d::coeffs<1,tl_real>, shapefcts2d::coeffs<1,tl_real> >
layer_potentials<1>( const std::array<point,3> X, geometry::tl_point P )
{
    using math::pifac;
    constexpr real eps { std::numeric_limits<real>::epsilon() }; 

       point E[3] = { X[1]-X[0], X[2]-X[1], X[0]-X[2] };
    tl_point R[3] = { X[0]-P,    X[1]-P,    X[2]-P    };
       point N    =  cross_prod(E[0],E[1]);
    const    real  e[3] = { E[0].r(),  E[1].r(),  E[2].r()  };
    const tl_real  r[3] = { R[0].r(), R[1].r(), R[2].r() };
    const    real  n    =   N.r();
    E[0] *= 1/e[0]; E[1] *= 1/e[1]; E[2] *= 1/e[2]; N *= 1/n;

    const tl_point RE[3] = { cross_prod(R[0],E[0]), cross_prod(R[1],E[1]), cross_prod(R[2],E[2]) };
    const tl_real  T     = scal_prod(R[0],cross_prod(R[1],R[2]));
    const tl_real  d     = scal_prod(R[0],N); // Signed distance of P to triangle's plane.

    const tl_real Omega = 2*atan2
    (
        T,
        r[0]*r[1]*r[2] +
        r[0]*scal_prod(R[1],R[2])+
        r[1]*scal_prod(R[0],R[2])+
        r[2]*scal_prod(R[0],R[1])
    );
  
    const tl_real L[3] =
    { 
         log( ( ( (r[1]+e[0])*(r[1]+e[0]) - r[0]*r[0]) + eps ) /
              ( (-(r[0]-e[0])*(r[0]-e[0]) + r[1]*r[1]) + eps ) ),

         log( ( ( (r[2]+e[1])*(r[2]+e[1]) - r[1]*r[1]) + eps ) /
              ( (-(r[1]-e[1])*(r[1]-e[1]) + r[2]*r[2]) + eps ) ),

         log( ( ( (r[0]+e[2])*(r[0]+e[2]) - r[2]*r[2]) + eps ) /
              ( (-(r[2]-e[2])*(r[2]-e[2]) + r[0]*r[0]) + eps ) )
    };

    const tl_real M[3]
    {
       (r[0]+r[1])*(e[0]+(r[1]-r[0])*(r[1]-r[0])/e[0])/4 + RE[0].r2()*L[0]/2,
       (r[1]+r[2])*(e[1]+(r[2]-r[1])*(r[2]-r[1])/e[1])/4 + RE[1].r2()*L[1]/2,
       (r[2]+r[0])*(e[2]+(r[0]-r[2])*(r[0]-r[2])/e[2])/4 + RE[2].r2()*L[2]/2
    };

    const tl_real  II  = scal_prod(N,RE[0]*L[0] + RE[1]*L[1] + RE[2]*L[2]) - d*Omega;
    const tl_point III = Omega*N + cross_prod(N,L[0]*E[0]+L[1]*E[1]+L[2]*E[2]);
    const tl_point IV  = II*N;
    const tl_point V   = M[0]*E[0] + M[1]*E[1] + M[2]*E[2];

    tl_point singlet = tl_point( scal_prod(RE[1], IV), scal_prod(RE[2], IV), scal_prod(RE[0] ,IV) ) -
                       tl_point( scal_prod( E[1],  V), scal_prod( E[2],  V), scal_prod( E[0],  V) );
    tl_point doublet = tl_point( scal_prod(RE[1],III), scal_prod(RE[2],III), scal_prod(RE[0],III) );
    singlet.x *= e[1]; singlet.y *= e[2]; singlet.z *= e[0];
    doublet.x *= e[1]; doublet.y *= e[2]; doublet.z *= e[0];
    singlet   *=  pifac/n;
    doublet   *= -pifac/n; 

    return std::pair< shapefcts2d::coeffs<1,tl_real>, shapefcts2d::coeffs<1,tl_real> >
    {
       {{ singlet.x, singlet.y, singlet.z }},
       {{ doublet.x, doublet.y, doublet.z }}
    };
}

template <>
point biot_savart<0>( const std::array<point,4> X, const std::array<point,1> OM, point P )
{  
    namespace bs = boost::simd;
    using bs::load;
    using bs::store;
    using pack_t = bs::pack<real>;

    using std::array;
    using math::pifac;

    constexpr size_t pack_size { bs::cardinal_of<pack_t>()            };
    constexpr real   eps       { std::numeric_limits<real>::epsilon() }; 
              real   work[16];

    //////////////////// TOPOLOGY INFORMATION //////////////////////
    // Face       Nodes                   Edges                   //
    // Face 0:    { X[1], X[2], X[3] }    {  E[2],  E[5], -E[4] } //
    // Face 1:    { X[0], X[3], X[2] }    {  E[3], -E[5], -E[1] } //
    // Face 2:    { X[0], X[1], X[3] }    {  E[0],  E[4], -E[3] } //
    // Face 3:    { X[0], X[2], X[1] }    {  E[1], -E[2], -E[0] } //
    ////////////////////////////////////////////////////////////////

    array<point,4> R  {{ X[0]-P, X[1]-P, X[2]-P, X[3]-P }};
    array<point,6> E  {{ X[1]-X[0], X[2]-X[0], X[2]-X[1],
                         X[3]-X[0], X[3]-X[1], X[3]-X[2] }};
    array<point,4> N  {{ cross_prod(E[2],E[5]), cross_prod(E[5],E[3]),
                         cross_prod(E[0],E[4]), cross_prod(E[2],E[1]) }};

    work[ 0] = E[0].r2();
    work[ 1] = E[1].r2();
    work[ 2] = E[2].r2();
    work[ 3] = E[3].r2();
    work[ 4] = E[4].r2();
    work[ 5] = E[5].r2();
    work[ 6] = R[0].r2();
    work[ 7] = R[1].r2();
    work[ 8] = R[2].r2();
    work[ 9] = R[3].r2();
    work[10] = N[0].r2();
    work[11] = N[1].r2();
    work[12] = N[2].r2();
    work[13] = N[3].r2();
    work[14] = 0;
    work[15] = 0;

    size_t i;
    for ( i = 0; i < 14 && i + pack_size <= 16; i += pack_size )
    {
        pack_t pack( work + i );
        bs::store( bs::sqrt(pack), work + i );
    }
    for ( ; i < 10; ++i ) work[i] = sqrt(work[i]);

    const array<real,6> e {{ work[ 0], work[ 1], work[ 2],
                             work[ 3], work[ 4], work[ 5] }};
    const array<real,4> r {{ work[ 6], work[ 7], work[ 8], work[ 9] }};
    const array<real,4> n {{ work[10], work[11], work[12], work[13] }};

    E[0] *= 1/e[0];
    E[1] *= 1/e[1];
    E[2] *= 1/e[2];
    E[3] *= 1/e[3];
    E[4] *= 1/e[4];
    E[5] *= 1/e[5];
    N[0] *= 1/n[0];
    N[1] *= 1/n[1];
    N[2] *= 1/n[2];
    N[3] *= 1/n[3];

    // Compute solid angles.
    work[0] = triple_prod(R[1],R[2],R[3]);
    work[1] = triple_prod(R[0],R[3],R[2]);
    work[2] = triple_prod(R[0],R[1],R[3]);
    work[3] = triple_prod(R[0],R[2],R[1]);
   
    work[4] = r[1]*r[2]*r[3] + r[1]*scal_prod(R[2],R[3]) + r[2]*scal_prod(R[1],R[3]) + r[3]*scal_prod(R[1],R[2]);
    work[5] = r[0]*r[3]*r[2] + r[0]*scal_prod(R[3],R[2]) + r[3]*scal_prod(R[0],R[2]) + r[2]*scal_prod(R[0],R[3]);
    work[6] = r[0]*r[1]*r[3] + r[0]*scal_prod(R[1],R[3]) + r[1]*scal_prod(R[0],R[3]) + r[3]*scal_prod(R[0],R[1]);
    work[7] = r[0]*r[2]*r[1] + r[0]*scal_prod(R[2],R[1]) + r[2]*scal_prod(R[0],R[1]) + r[1]*scal_prod(R[0],R[2]);
    for ( i = 0; i + pack_size <= 4; i += pack_size )
    {
        pack_t  y(work + i);
        pack_t  x(work + i + 4 );
        bs::store( bs::pedantic_(bs::atan2)(y,x), work + i );
    }
    for ( ; i < 4; ++i ) work[i] = atan2(work[i],work[i+4]);
    const array<real,4> Omega {{ 2*work[0], 2*work[1], 2*work[2], 2*work[3] }};


    // Compute L terms.
    work[0 ] = (r[1]+e[0])*(r[1]+e[0]) - r[0]*r[0];
    work[1 ] = (r[2]+e[1])*(r[2]+e[1]) - r[0]*r[0];
    work[2 ] = (r[2]+e[2])*(r[2]+e[2]) - r[1]*r[1];
    work[3 ] = (r[3]+e[3])*(r[3]+e[3]) - r[0]*r[0];
    work[4 ] = (r[3]+e[4])*(r[3]+e[4]) - r[1]*r[1];
    work[5 ] = (r[3]+e[5])*(r[3]+e[5]) - r[2]*r[2];
    work[6 ] = 1;
    work[7 ] = 1;

    work[8 ] = -(r[0]-e[0])*(r[0]-e[0]) + r[1]*r[1];
    work[9 ] = -(r[0]-e[1])*(r[0]-e[1]) + r[2]*r[2];
    work[10] = -(r[1]-e[2])*(r[1]-e[2]) + r[2]*r[2];
    work[11] = -(r[0]-e[3])*(r[0]-e[3]) + r[3]*r[3];
    work[12] = -(r[1]-e[4])*(r[1]-e[4]) + r[3]*r[3];
    work[13] = -(r[2]-e[5])*(r[2]-e[5]) + r[3]*r[3];
    work[14] = 1;
    work[15] = 1;

    for ( i = 0; i < 6 && i + pack_size <= 8; i += pack_size )
    {
        pack_t   numerator( work + i );
        pack_t denominator( work + i + 8 );
        numerator += eps; denominator += eps;
        pack_t args = numerator / denominator;
        bs::store( bs::log(args), work + i );
    }
    for ( ; i < 6; ++i ) work[i] = log((work[i]+eps)/(work[i+8]+eps));
    const array<real,6> L {{ work[0], work[1], work[2], work[3], work[4], work[5] }};
  

    #define COMPUTE_SIGMA( NF, R0, R1, R2, E0, E1, E2, L0, L1, L2, OmegaF, result )     \
    {                                                                                   \
        const point RE[3] = { cross_prod((R0),(E0)), cross_prod((R1),(E1)),             \
                              cross_prod((R2),(E2)) };                                  \
        const real d = scal_prod((R0),(NF));                                            \
        result = scal_prod((NF), RE[0]*(L0) + RE[1]*(L1) + RE[2]*(L2)) - (d)*(OmegaF);  \
    }
    
    real sigma[4];
    COMPUTE_SIGMA( N[0], R[1], R[2], R[3], E[2], E[5],-E[4], L[2], L[5], L[4], Omega[0], sigma[0] );
    COMPUTE_SIGMA( N[1], R[0], R[3], R[2], E[3],-E[5],-E[1], L[3], L[5], L[1], Omega[1], sigma[1] );
    COMPUTE_SIGMA( N[2], R[0], R[1], R[3], E[0], E[4],-E[3], L[0], L[4], L[3], Omega[2], sigma[2] );
    COMPUTE_SIGMA( N[3], R[0], R[2], R[1], E[1],-E[2],-E[0], L[1], L[2], L[0], Omega[3], sigma[3] );
    #undef COMPUTE_SIGMA

    point result { cross_prod(OM[0],N[0])*sigma[0] + cross_prod(OM[0],N[1])*sigma[1] +
                   cross_prod(OM[0],N[2])*sigma[2] + cross_prod(OM[0],N[3])*sigma[3] };
    result *= pifac;
        
    const bool negative { scal_prod(E[0],N[0]) < 0 };
    return negative ? -result : result;
}

template <>
point biot_savart<1>( const std::array<point,4> X, const std::array<point,4> OM, point P )
{  
    namespace bs = boost::simd;
    using bs::load;
    using bs::store;
    using pack_t = bs::pack<real>;

    using std::array;
    using math::pifac;

    constexpr size_t pack_size { bs::cardinal_of<pack_t>()            };
    constexpr real   eps       { std::numeric_limits<real>::epsilon() }; 
              real   work[16];

    const tensor grad_omega = ( tensor { (OM[1]-OM[0]).x, (OM[2]-OM[0]).x, (OM[3]-OM[0]).x,
                                         (OM[1]-OM[0]).y, (OM[2]-OM[0]).y, (OM[3]-OM[0]).y,
                                         (OM[1]-OM[0]).z, (OM[2]-OM[0]).z, (OM[3]-OM[0]).z } ) *
                              ( tensor { (X [1]-X [0]).x, (X [2]-X [0]).x, (X[3]-X[0]).x,
                                         (X [1]-X [0]).y, (X [2]-X [0]).y, (X[3]-X[0]).y,
                                         (X [1]-X [0]).z, (X [2]-X [0]).z, (X[3]-X[0]).z } ).inv();
    const point halfrot_omega = point( grad_omega(2,1) - grad_omega(1,2),
                                       grad_omega(0,2) - grad_omega(2,0),
                                       grad_omega(1,0) - grad_omega(0,1) )/2;

    //////////////////// TOPOLOGY INFORMATION //////////////////////
    // Face       Nodes                   Edges                   //
    // Face 0:    { X[1], X[2], X[3] }    {  E[2],  E[5], -E[4] } //
    // Face 1:    { X[0], X[3], X[2] }    {  E[3], -E[5], -E[1] } //
    // Face 2:    { X[0], X[1], X[3] }    {  E[0],  E[4], -E[3] } //
    // Face 3:    { X[0], X[2], X[1] }    {  E[1], -E[2], -E[0] } //
    ////////////////////////////////////////////////////////////////
    array<point,4> R  {{ X[0]-P, X[1]-P, X[2]-P, X[3]-P }};
    array<point,6> E  {{ X[1]-X[0], X[2]-X[0], X[2]-X[1],
                         X[3]-X[0], X[3]-X[1], X[3]-X[2] }};
    array<point,4> N  {{ cross_prod(E[2],E[5]), cross_prod(E[5],E[3]),
                         cross_prod(E[0],E[4]), cross_prod(E[2],E[1]) }};

    work[ 0] = E[0].r2();
    work[ 1] = E[1].r2();
    work[ 2] = E[2].r2();
    work[ 3] = E[3].r2();
    work[ 4] = E[4].r2();
    work[ 5] = E[5].r2();
    work[ 6] = R[0].r2();
    work[ 7] = R[1].r2();
    work[ 8] = R[2].r2();
    work[ 9] = R[3].r2();
    work[10] = N[0].r2();
    work[11] = N[1].r2();
    work[12] = N[2].r2();
    work[13] = N[3].r2();
    work[14] = 0;
    work[15] = 0;

    size_t i;
    for ( i = 0; i < 14 && i + pack_size <= 16; i += pack_size )
    {
        pack_t pack( work + i );
        bs::store( bs::sqrt(pack), work + i );
    }
    for ( ; i < 14; ++i ) work[i] = sqrt(work[i]);

    const array<real,6> e {{ work[ 0], work[ 1], work[ 2],
                             work[ 3], work[ 4], work[ 5] }};
    const array<real,4> r {{ work[ 6], work[ 7], work[ 8], work[ 9] }};
    const array<real,4> n {{ work[10], work[11], work[12], work[13] }};

    E[0] *= 1/e[0];
    E[1] *= 1/e[1];
    E[2] *= 1/e[2];
    E[3] *= 1/e[3];
    E[4] *= 1/e[4];
    E[5] *= 1/e[5];
    N[0] *= 1/n[0];
    N[1] *= 1/n[1];
    N[2] *= 1/n[2];
    N[3] *= 1/n[3];

    // Compute solid angles.
    work[0] = triple_prod(R[1],R[2],R[3]);
    work[1] = triple_prod(R[0],R[3],R[2]);
    work[2] = triple_prod(R[0],R[1],R[3]);
    work[3] = triple_prod(R[0],R[2],R[1]);
   
    work[4] = r[1]*r[2]*r[3] + r[1]*scal_prod(R[2],R[3]) + r[2]*scal_prod(R[1],R[3]) + r[3]*scal_prod(R[1],R[2]);
    work[5] = r[0]*r[3]*r[2] + r[0]*scal_prod(R[3],R[2]) + r[3]*scal_prod(R[0],R[2]) + r[2]*scal_prod(R[0],R[3]);
    work[6] = r[0]*r[1]*r[3] + r[0]*scal_prod(R[1],R[3]) + r[1]*scal_prod(R[0],R[3]) + r[3]*scal_prod(R[0],R[1]);
    work[7] = r[0]*r[2]*r[1] + r[0]*scal_prod(R[2],R[1]) + r[2]*scal_prod(R[0],R[1]) + r[1]*scal_prod(R[0],R[2]);
    for ( i = 0; i + pack_size <= 4; i += pack_size )
    {
        pack_t  y(work + i);
        pack_t  x(work + i + 4 );
        bs::store( bs::pedantic_(bs::atan2)(y,x), work + i );
    }
    for ( ; i < 4; ++i ) work[i] = atan2(work[i],work[i+4]);
    const array<real,4> Omega {{ 2*work[0], 2*work[1], 2*work[2], 2*work[3] }};

    // Compute L and M terms.
    work[0 ] = (r[1]+e[0])*(r[1]+e[0]) - r[0]*r[0];
    work[1 ] = (r[2]+e[1])*(r[2]+e[1]) - r[0]*r[0];
    work[2 ] = (r[2]+e[2])*(r[2]+e[2]) - r[1]*r[1];
    work[3 ] = (r[3]+e[3])*(r[3]+e[3]) - r[0]*r[0];
    work[4 ] = (r[3]+e[4])*(r[3]+e[4]) - r[1]*r[1];
    work[5 ] = (r[3]+e[5])*(r[3]+e[5]) - r[2]*r[2];
    work[6 ] = 1;
    work[7 ] = 1;

    work[8 ] = -(r[0]-e[0])*(r[0]-e[0]) + r[1]*r[1];
    work[9 ] = -(r[0]-e[1])*(r[0]-e[1]) + r[2]*r[2];
    work[10] = -(r[1]-e[2])*(r[1]-e[2]) + r[2]*r[2];
    work[11] = -(r[0]-e[3])*(r[0]-e[3]) + r[3]*r[3];
    work[12] = -(r[1]-e[4])*(r[1]-e[4]) + r[3]*r[3];
    work[13] = -(r[2]-e[5])*(r[2]-e[5]) + r[3]*r[3];
    work[14] = 1;
    work[15] = 1;

    for ( i = 0; i < 6 && i + pack_size <= 16; i += pack_size )
    {
        pack_t   numerator( work + i );
        pack_t denominator( work + i + 8 );
        numerator += eps; denominator += eps;
        pack_t args = numerator / denominator;
        bs::store( bs::log(args), work + i );
    }
    for ( ; i < 6; ++i ) work[i] = log((work[i]+eps)/(work[i+8]+eps));
    const array<real,6> L {{ work[0], work[1], work[2], work[3], work[4], work[5] }};
    const array<real,6> M
    {{
        (r[1]+r[0])*(e[0]+(r[1]-r[0])*(r[1]-r[0])/e[0])/4 + cross_prod(R[1],E[0]).r2()*L[0]/2,
        (r[2]+r[0])*(e[1]+(r[2]-r[0])*(r[2]-r[0])/e[1])/4 + cross_prod(R[2],E[1]).r2()*L[1]/2,
        (r[2]+r[1])*(e[2]+(r[2]-r[1])*(r[2]-r[1])/e[2])/4 + cross_prod(R[2],E[2]).r2()*L[2]/2,
        (r[3]+r[0])*(e[3]+(r[3]-r[0])*(r[3]-r[0])/e[3])/4 + cross_prod(R[3],E[3]).r2()*L[3]/2,
        (r[3]+r[1])*(e[4]+(r[3]-r[1])*(r[3]-r[1])/e[4])/4 + cross_prod(R[3],E[4]).r2()*L[4]/2,
        (r[3]+r[2])*(e[5]+(r[3]-r[2])*(r[3]-r[2])/e[5])/4 + cross_prod(R[3],E[5]).r2()*L[5]/2
    }};

    #define COMPUTE_POT( R0, R1, R2, E0, E1, E2, e0, e1, e2, L0, L1, L2, M0, M1, M2, NF, nf, OmegaF, singlet ) \
    {                                                                                              \
        const point RE[3] = { cross_prod((R0),(E0)), cross_prod((R1),(E1)),                        \
                              cross_prod((R2),(E2)) };                                             \
        const real d = scal_prod((R0),(NF));                                                       \
        const real  II = scal_prod((NF), RE[0]*(L0) + RE[1]*(L1) + RE[2]*(L2)) - (d)*(OmegaF);     \
        const point IV = II*NF;                                                                    \
        const point  V = (E0)*(M0) + (E1)*(M1) + (E2)*(M2);                                        \
        singlet = point( scal_prod(RE[1], IV), scal_prod(RE[2], IV), scal_prod(RE[0], IV) ) -      \
                  point( scal_prod( (E1),  V), scal_prod( (E2),  V), scal_prod( (E0),  V) );       \
        singlet.x *= e1; singlet.y *= e2; singlet.z *= e0;                                         \
        singlet /= nf; \
    }
    

    // Face 0.
    point POT, velocity, sigma0, sigma1, sigma2;
    COMPUTE_POT(   R[1] ,R[2], R[3], 
                   E[2], E[5],-E[4], e[2], e[5], e[4],
                   L[2], L[5], L[4], M[2], M[5], M[4],
                   N[0], n[0], Omega[0], POT );
    sigma0 = point { halfrot_omega*scal_prod(R[1],N[0]) + cross_prod(OM[1],N[0]) };
    sigma1 = point { halfrot_omega*scal_prod(R[2],N[0]) + cross_prod(OM[2],N[0]) };
    sigma2 = point { halfrot_omega*scal_prod(R[3],N[0]) + cross_prod(OM[3],N[0]) };
    velocity += sigma0*POT.x + sigma1*POT.y + sigma2*POT.z;

    // Face 1.
    COMPUTE_POT(   R[0], R[3], R[2], 
                   E[3],-E[5],-E[1], e[3], e[5], e[1],
                   L[3], L[5], L[1], M[3], M[5], M[1],
                   N[1], n[1], Omega[1], POT );
    sigma0 = point { halfrot_omega*scal_prod(R[0],N[1]) + cross_prod(OM[0],N[1]) };
    sigma1 = point { halfrot_omega*scal_prod(R[3],N[1]) + cross_prod(OM[3],N[1]) };
    sigma2 = point { halfrot_omega*scal_prod(R[2],N[1]) + cross_prod(OM[2],N[1]) };
    velocity += sigma0*POT.x + sigma1*POT.y + sigma2*POT.z;

    // Face 2.
    COMPUTE_POT(   R[0], R[1], R[3], 
                   E[0], E[4],-E[3], e[0], e[4], e[3],
                   L[0], L[4], L[3], M[0], M[4], M[3],
                   N[2], n[2], Omega[2], POT );
    sigma0 = point { halfrot_omega*scal_prod(R[0],N[2]) + cross_prod(OM[0],N[2]) };
    sigma1 = point { halfrot_omega*scal_prod(R[1],N[2]) + cross_prod(OM[1],N[2]) };
    sigma2 = point { halfrot_omega*scal_prod(R[3],N[2]) + cross_prod(OM[3],N[2]) };
    velocity += sigma0*POT.x + sigma1*POT.y + sigma2*POT.z;

    // Face 3.
    COMPUTE_POT(   R[0], R[2], R[1], 
                   E[1],-E[2],-E[0], e[1], e[2], e[0],
                   L[1], L[2], L[0], M[1], M[2], M[0],
                   N[3], n[3], Omega[3], POT );
    sigma0 = point { halfrot_omega*scal_prod(R[0],N[3]) + cross_prod(OM[0],N[3]) };
    sigma1 = point { halfrot_omega*scal_prod(R[2],N[3]) + cross_prod(OM[2],N[3]) };
    sigma2 = point { halfrot_omega*scal_prod(R[1],N[3]) + cross_prod(OM[1],N[3]) };
    velocity += sigma0*POT.x + sigma1*POT.y + sigma2*POT.z;
    #undef COMPUTE_POT

    velocity *= pifac;
    const bool negative { scal_prod(E[0],N[0]) < 0 };
    return  negative ? -velocity : velocity;
}

}

