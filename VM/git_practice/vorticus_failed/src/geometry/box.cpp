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
#include "geometry/box.h"

namespace geometry
{

/*!
 * \brief Tests if a triangle and a box intersect.
 * \see    Michael Schwarz and Hans-Peter Seidelm
 *         "Fast Parallel Surface and Solid Voxelization on GPUs", 
 *         ACM Transactions on Graphics, 2010.
 *
 * This is an implementation of the test by Schwarz and Seidel.
 */
bool intersects_triangle( const box &B, const std::array<point,3> &T ) noexcept
{
    using std::max;

    box BT { T[0], T[0] };
    BT = bounding_box(BT,T[1]);
    BT = bounding_box(BT,T[2]);
    

    point  p = B.min;
    point dp = B.max - B.min;
    point e[3] { T[1]-T[0], T[2]-T[1], T[0]-T[2] };

    point n = cross_prod( e[0], e[1] );
    point c { (n.x > 0) ? dp.x : 0,
              (n.y > 0) ? dp.y : 0,
              (n.z > 0) ? dp.z : 0 };
    real d1 = scal_prod( n,     c - T[0] );
    real d2 = scal_prod( n, (dp-c)- T[0] );

    // xy.
    real ne_xy[3][2];
    if ( n.z < 0 )
    {
        ne_xy[0][0] =  e[0].y; ne_xy[0][1] = -e[0].x;
        ne_xy[1][0] =  e[1].y; ne_xy[1][1] = -e[1].x;
        ne_xy[2][0] =  e[2].y; ne_xy[2][1] = -e[2].x;
    }
    else
    {
        ne_xy[0][0] = -e[0].y; ne_xy[0][1] =  e[0].x;
        ne_xy[1][0] = -e[1].y; ne_xy[1][1] =  e[1].x;
        ne_xy[2][0] = -e[2].y; ne_xy[2][1] =  e[2].x;
    }

    // xz.
    real ne_xz[3][2];
    if ( n.y > 0 )
    {
        ne_xz[0][0] =  e[0].z; ne_xz[0][1] = -e[0].x;
        ne_xz[1][0] =  e[1].z; ne_xz[1][1] = -e[1].x;
        ne_xz[2][0] =  e[2].z; ne_xz[2][1] = -e[2].x;
    }
    else
    {
        ne_xz[0][0] = -e[0].z; ne_xz[0][1] =  e[0].x;
        ne_xz[1][0] = -e[1].z; ne_xz[1][1] =  e[1].x;
        ne_xz[2][0] = -e[2].z; ne_xz[2][1] =  e[2].x;
    }

    // yz.
    real ne_yz[3][2];
    if ( n.x < 0 )
    {
        ne_yz[0][0] =  e[0].z; ne_yz[0][1] = -e[0].y;
        ne_yz[1][0] =  e[1].z; ne_yz[1][1] = -e[1].y;
        ne_yz[2][0] =  e[2].z; ne_yz[2][1] = -e[2].y;
    }
    else
    {
        ne_yz[0][0] = -e[0].z; ne_yz[0][1] =  e[0].y;
        ne_yz[1][0] = -e[1].z; ne_yz[1][1] =  e[1].y;
        ne_yz[2][0] = -e[2].z; ne_yz[2][1] =  e[2].y;
    }

    real de_xy[3];
    de_xy[0] = -(ne_xy[0][0]*T[0].x + ne_xy[0][1]*T[0].y) + max(0.,ne_xy[0][0]*dp.x) + max(0.,ne_xy[0][1]*dp.y);
    de_xy[1] = -(ne_xy[1][0]*T[1].x + ne_xy[1][1]*T[1].y) + max(0.,ne_xy[1][0]*dp.x) + max(0.,ne_xy[1][1]*dp.y);
    de_xy[2] = -(ne_xy[2][0]*T[2].x + ne_xy[2][1]*T[2].y) + max(0.,ne_xy[2][0]*dp.x) + max(0.,ne_xy[2][1]*dp.y);

    real de_xz[3];
    de_xz[0] = -(ne_xz[0][0]*T[0].x + ne_xz[0][1]*T[0].z) + max(0.,ne_xz[0][0]*dp.x) + max(0.,ne_xz[0][1]*dp.z);
    de_xz[1] = -(ne_xz[1][0]*T[1].x + ne_xz[1][1]*T[1].z) + max(0.,ne_xz[1][0]*dp.x) + max(0.,ne_xz[1][1]*dp.z);
    de_xz[2] = -(ne_xz[2][0]*T[2].x + ne_xz[2][1]*T[2].z) + max(0.,ne_xz[2][0]*dp.x) + max(0.,ne_xz[2][1]*dp.z);

    real de_yz[3];
    de_yz[0] = -(ne_yz[0][0]*T[0].y + ne_yz[0][1]*T[0].z) + max(0.,ne_yz[0][0]*dp.y) + max(0.,ne_yz[0][1]*dp.z);
    de_yz[1] = -(ne_yz[1][0]*T[1].y + ne_yz[1][1]*T[1].z) + max(0.,ne_yz[1][0]*dp.y) + max(0.,ne_yz[1][1]*dp.z);
    de_yz[2] = -(ne_yz[2][0]*T[2].y + ne_yz[2][1]*T[2].z) + max(0.,ne_yz[2][0]*dp.y) + max(0.,ne_yz[2][1]*dp.z);

    real np { scal_prod(n,p) };
    bool test_box   = (dist(B,BT) == 0); 
    bool test_plane = ( ((np + d1)*(np + d2)) <= 0 );
    bool test_xy  = ( (ne_xy[0][0]*p.x + ne_xy[0][1]*p.y + de_xy[0]) >= 0 ) &&
                    ( (ne_xy[1][0]*p.x + ne_xy[1][1]*p.y + de_xy[1]) >= 0 ) &&
                    ( (ne_xy[2][0]*p.x + ne_xy[2][1]*p.y + de_xy[2]) >= 0 );

    bool test_xz  = ( (ne_xz[0][0]*p.x + ne_xz[0][1]*p.z + de_xz[0]) >= 0 ) &&
                    ( (ne_xz[1][0]*p.x + ne_xz[1][1]*p.z + de_xz[1]) >= 0 ) &&
                    ( (ne_xz[2][0]*p.x + ne_xz[2][1]*p.z + de_xz[2]) >= 0 );

    bool test_yz  = ( (ne_yz[0][0]*p.y + ne_yz[0][1]*p.z + de_yz[0]) >= 0 ) &&
                    ( (ne_yz[1][0]*p.y + ne_yz[1][1]*p.z + de_yz[1]) >= 0 ) &&
                    ( (ne_yz[2][0]*p.y + ne_yz[2][1]*p.z + de_yz[2]) >= 0 );

    return test_box && test_plane && test_xy && test_xz && test_yz;
}

namespace
{

constexpr
point support( const box &B, const point &d )
{
    return point
    {
        ( d.x < 0 ) ? B.min.x : B.max.x,
        ( d.y < 0 ) ? B.min.y : B.max.y,
        ( d.z < 0 ) ? B.min.z : B.max.z
    };
}

template <size_t N> inline
point support( const std::array<point,N> &P, const point &d )
{
    size_t n_max { 0 };
    real   e_max { scal_prod(P[0],d) };
    for ( size_t n = 1; n < N; ++n )
    {
        real e = scal_prod(P[n],d);
        if ( e > e_max )
        {
            e_max = e;
            n_max = n;
        }
    }
    return P[n_max];
}

template <size_t NN>
void do_simplex( std::array<point,4> &S, size_t &N, point &d );

template <> inline
void do_simplex<2>( std::array<point,4> &S, size_t &N, point &d )
{
    const point A  { S[1] };
    const point B  { S[0] };
    const point AB { B-A  };
    #define GJK_TEST(x) (scal_prod((x),A) < 0)
    if ( GJK_TEST(B-A) )
    {
        N = 2;
        d = cross_prod( cross_prod(A,AB), AB);
    }
    else
    {
        N = 1;
        S[0] = A;
        d = -A;
    }
    #undef  GJK_TEST
}

template <> inline
void do_simplex<3>( std::array<point,4> &S, size_t &N, point &d )
{
    const point A { S[2] };
    const point B { S[1] };
    const point C { S[0] };
    const point AB { B-A };
    const point AC { C-A };
    const point ABC { cross_prod(AB,AC) };

    #define GJK_TEST(x) (scal_prod((x),A) < 0)
    if ( GJK_TEST(cross_prod(ABC,AC)) )
    {
       if ( GJK_TEST(AC) )
       {
            N = 2; // S = { A, C }
            S[1] = A;
            S[0] = C;
            d = cross_prod(cross_prod(S[0],AC),AC);
       }
       else
       {
            N = 2; // S = { A, B }
            S[1] = A;
            S[0] = B;
            do_simplex<2>(S,N,d);
       }
    } 
    else if ( GJK_TEST(cross_prod(AB,ABC)) )
    {
            N = 2; // S = { A, B }
            S[1] = A;
            S[0] = B;
            do_simplex<2>(S,N,d);
    }
    else if ( GJK_TEST(ABC) )
    {
        N = 3; // S = { A, B, C }
        d = ABC;
    }
    else
    {
        N = 3; // S = { A, C, B }
        S[2] = A;
        S[1] = C;
        S[0] = B;
        d = -ABC;
    }
    #undef GJK_TEST
}

template <> inline
void do_simplex<4>( std::array<point,4> &S, size_t &N, point &d )
{
    const point A { S[3] };
    const point B { S[2] };
    const point C { S[1] };
    const point D { S[0] };
    const point AB { B-A };
    const point AC { C-A };
    const point AD { D-A };
 
    #define GJK_TEST(x) (scal_prod((x),A) < 0)
    if ( GJK_TEST(cross_prod(AB,AC)) )
    {
        N = 3; // S = { A, B, C };
        S[2] = A;
        S[1] = B;
        S[0] = C;
    }
    else if ( GJK_TEST(cross_prod(AC,AD)) )
    {
        N = 3; // S = { A, C, D };
        S[2] = A;
        S[1] = C;
        S[0] = D;
    }
    else if ( GJK_TEST(cross_prod(AD,AB)) )
    {
        N = 3; // S = { A, D, B };
        S[2] = A;
        S[1] = D;
        S[0] = B;
    }
    else return;

    do_simplex<3>(S,N,d);
    #undef GJK_TEST
}

}

// Detect intersection using the Gilbert-Johnson-Keerthi (GJK) algorithm.
// See Gino van den Bergen, "Collision Detection in Interactive 3D Environments",
// Algorithm 4.5.
bool intersects_tetrahedron( const box &B, const std::array<point,4> &T ) noexcept
{
    constexpr real eps_tol { 128.*std::numeric_limits<real>::epsilon() };

    point v {1,0,0};
    std::array<point,4> W; size_t NW {0}; real w_max {0};
    std::array<point,4> Y; size_t NY {0};

    while ( NW != 4 && (v.r2() > eps_tol*w_max) )
    {
        point w = support(B,-v) - support(T,v);

        switch (NY)
        {
        case 4:  if ( w == Y[3] ) return true; // Fall-through!
        case 3:  if ( w == Y[2] ) return true;
        case 2:  if ( w == Y[1] ) return true;
        case 1:  if ( w == Y[0] ) return true;
        default: if ( scal_prod(w,v) > 0 ) return false;
        }

        Y = W; NY = NW;
        Y[ NY++ ] = w;
        W = Y; NW = NY;

        switch (NW)
        {
        case 4: do_simplex<4>( W, NW, v ); break;
        case 3: do_simplex<3>( W, NW, v ); break;
        case 2: do_simplex<2>( W, NW, v ); break;
        }

        w_max = 0;
        switch (NW)
        {
        case 4: w_max = std::max( w_max, W[3].r2() ); // Fall-through!
        case 3: w_max = std::max( w_max, W[2].r2() ); 
        case 2: w_max = std::max( w_max, W[1].r2() );
        case 1: w_max = std::max( w_max, W[0].r2() );
        }
    }
    return true;
}

}

