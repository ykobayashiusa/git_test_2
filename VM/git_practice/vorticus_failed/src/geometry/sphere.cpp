/*
 * Copyright (C) 2017 Matthias Kirchhart
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
#include "geometry/sphere.h" 

namespace geometry
{

namespace
{

/*!
 * \brief Solve systems arising in method mb in the least-squares sense.
 * 
 * In the method "mb()" one needs to apply the Moore-Penrose pseudoinverse
 * of the resulting matrix Q^T to two right-hand sides. This method achieves 
 * this. Expects A = Q^T as a C-style array. Assumes that A has full rank.
 */
template <size_t n>
void solve( const real A[n][3], const real a[n], const real b[n], point &e, point &f );

template <> inline
void solve<1>( const real A[1][3], const real a[1], const real b[1],
               point &e, point &f )
{
    const real inv_sigma { 1./(A[0][0]*A[0][0] + A[0][1]*A[0][1] + A[0][2]*A[0][2]) };
    const real fac_a { a[0]*inv_sigma };
    const real fac_b { b[0]*inv_sigma };
    e.x = A[0][0]*fac_a; f.x = A[0][0]*fac_b;
    e.y = A[0][1]*fac_a; f.y = A[0][1]*fac_b;
    e.z = A[0][2]*fac_a; f.z = A[0][2]*fac_b;
}

template <> inline
void solve<2>( const real A[2][3], const real a[2], const real b[2],
               point &e, point &f )
{
    // Matrix A*A^T.
    const real AAt[2][2]
    {
        { A[0][0]*A[0][0]+A[0][1]*A[0][1]+A[0][2]*A[0][2], A[0][0]*A[1][0]+A[0][1]*A[1][1]+A[0][2]*A[1][2] },
        { A[0][0]*A[1][0]+A[0][1]*A[1][1]+A[0][2]*A[1][2], A[1][0]*A[1][0]+A[1][1]*A[1][1]+A[1][2]*A[1][2] }
    };
    const real det_inv { 1./(AAt[0][0]*AAt[1][1] - AAt[0][1]*AAt[1][0]) };

    const real a_tmp[2]
    {
        det_inv*(AAt[1][1]*a[0] - AAt[0][1]*a[1]),
        det_inv*(AAt[0][0]*a[1] - AAt[1][0]*a[0])
    };

    const real b_tmp[2]
    {
        det_inv*(AAt[1][1]*b[0] - AAt[0][1]*b[1]),
        det_inv*(AAt[0][0]*b[1] - AAt[1][0]*b[0])
    };

    e.x = A[0][0]*a_tmp[0] + A[1][0]*a_tmp[1];
    e.y = A[0][1]*a_tmp[0] + A[1][1]*a_tmp[1];
    e.z = A[0][2]*a_tmp[0] + A[1][2]*a_tmp[1];

    f.x = A[0][0]*b_tmp[0] + A[1][0]*b_tmp[1];
    f.y = A[0][1]*b_tmp[0] + A[1][1]*b_tmp[1];
    f.z = A[0][2]*b_tmp[0] + A[1][2]*b_tmp[1];
}

template <> inline
void solve<3>( const real A[3][3], const real a[3], const real b[3],
               point &e, point &f )
{
    // This is actually the easiest case, as we can take the actual
    // inverse of A. First compute the adjoint matrix.
    const real adj[3][3]
    {
        { A[1][1]*A[2][2] - A[1][2]*A[2][1], A[0][2]*A[2][1] - A[0][1]*A[2][2], A[0][1]*A[1][2] - A[0][2]*A[1][1] },
        { A[1][2]*A[2][0] - A[1][0]*A[2][2], A[0][0]*A[2][2] - A[0][2]*A[2][0], A[0][2]*A[1][0] - A[0][0]*A[1][2] },
        { A[1][0]*A[2][1] - A[1][1]*A[2][0], A[0][1]*A[2][0] - A[0][0]*A[2][1], A[0][0]*A[1][1] - A[0][1]*A[1][0] }
    };
    const real det_inv = 1./
    ( 
        A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] -
        A[2][0]*A[1][1]*A[0][2] - A[2][1]*A[1][2]*A[0][0] - A[2][2]*A[1][0]*A[0][1]
    );

    e.x = det_inv*(adj[0][0]*a[0] + adj[0][1]*a[1] + adj[0][2]*a[2]);
    e.y = det_inv*(adj[1][0]*a[0] + adj[1][1]*a[1] + adj[1][2]*a[2]);
    e.z = det_inv*(adj[2][0]*a[0] + adj[2][1]*a[1] + adj[2][2]*a[2]);
    
    f.x = det_inv*(adj[0][0]*b[0] + adj[0][1]*b[1] + adj[0][2]*b[2]);
    f.y = det_inv*(adj[1][0]*b[0] + adj[1][1]*b[1] + adj[1][2]*b[2]);
    f.z = det_inv*(adj[2][0]*b[0] + adj[2][1]*b[1] + adj[2][2]*b[2]);
}

/*!
 * \brief Given a basis V, combute its corresponding bounding sphere.
 * \param V The basis.
 * \return The resulting sphere.
 * \see Kaspar Fischer and Bernd Gärtner, "The Smallest Enclosing Ball of Balls:
 *      Combinatorial Structure and Algorithms", Lemma 3.5.
 *
 * Given a set of n + 1 spheres in three-dimensional space, n <= 3, this method
 * computes the minimal bounding sphere. The function assumes that the given
 * spheres form a basis, i.e., the bounding sphere of a subset of V does not
 * already contain the remaining spheres. Otherwise the result will be garbage.
 */
template <size_t n> inline
sphere mb( const std::array<sphere,n+1> &V )
{
    static_assert ( 0 < n && n < 4, "Invalid n in basis-computation." );

    // Matrix Q^T as a C-style matrix (row-major).
    real Qt[n][3], a[n], b[n], max_rho { V[n].radius };
    const real rn { V[n].radius };
    for ( size_t i = 0; i < n; ++i )
    {
        Qt[i][0] = V[i].centre.x - V[n].centre.x;
        Qt[i][1] = V[i].centre.y - V[n].centre.y;
        Qt[i][2] = V[i].centre.z - V[n].centre.z;
        
        const real ri { V[i].radius };
        a[i] = ri - rn;
        b[i] = ((V[i].centre-V[n].centre).r2() + rn*rn - ri*ri)/2;
        max_rho = std::max( max_rho, ri );
    }
    point e, f; solve<n>( Qt, a, b, e, f );

    // At this point we know the centre of the desired circle is
    // V[n].centre + f + rho*e. It remains to find rho.
    const real P = (rn + scal_prod(e,f))/(1-e.r2());
    const real Q = (f.r2() - rn*rn)/(1-e.r2());
    const real root = std::sqrt(std::abs(P*P+Q));
    const real rho  = ( (P-root) >= max_rho ) ? (P-root) : (P+root);

    return sphere { V[n].centre + rho*e + f, rho };
}

template <> inline
sphere mb<0>( const std::array<sphere,1> &V )
{
    return V[0];
}


/*!
 * \brief Check if a candidate basis VC does not violate the old given basis V.
 * \param V  The old basis.
 * \param n  The size of the old basis.
 * \param VB The ball corresponding to the old basis.
 * \param VC The new candidate basis.
 * \param nc The size of the new candidate basis.
 * \return Returns true if the the ball spanned by VC contains all balls of
 *         the old basis V, false otherwise.
 *
 * We are given a basis V and its corresponding bounding ball VB. Given a new
 * candidate basis VC, this method checks if the ball spanned by VC contains
 * all balls of the old basis. If so, the old basis and its corresponding
 * ball are set to the candidate basis and its corresponding ball and the
 * method returns true. Otherwise V and its ball remain unchanged and the method
 * returns false.
 */
bool try_basis(       std::array<sphere,4> &V,         size_t &n, sphere &VB,
                const std::array<sphere,4> &VC, const  size_t  nc )
{
    sphere B_candidate;
    switch ( nc )
    {
    case 1: B_candidate = mb<0>( {{ VC[0]                      }} ); break;
    case 2: B_candidate = mb<1>( {{ VC[0], VC[1]               }} ); break;
    case 3: B_candidate = mb<2>( {{ VC[0], VC[1], VC[2]        }} ); break;
    case 4: B_candidate = mb<3>( {{ VC[0], VC[1], VC[2], VC[3] }} ); break;
    }

    for ( size_t i = 0; i < n; ++i )
    {
        if ( ! B_candidate.approx_contains(V[i]) )
        {
            return false;
        }
    }

    V = VC; n = nc; VB = B_candidate;
    return true;
}

void basis_computation_1( std::array<sphere,4> &V, size_t &n, sphere &VB, sphere B )
{
    // n = 1.
    constexpr size_t N { 2 };
    const std::array<sphere,4> VCs[N] {
        { B },
        { V[0], B }
    };
    const size_t NCs[N] { 1, 2 };

    for ( size_t i = 0; i < N; ++i )
    {
        if ( try_basis(V,n,VB,VCs[i],NCs[i]) )
            return;
    }
    throw std::runtime_error { "Error while computing bounding sphere. No basis found. Case n = 1." };
}

void basis_computation_2( std::array<sphere,4> &V, size_t &n, sphere &VB, sphere B )
{
    // n = 2.
    constexpr size_t N {4};
    const std::array<sphere,4> VCs[N]
    {
        { B },
        { V[0], B },
        { V[1], B },
        { V[0], V[1], B }
    };
    const size_t NCs[N] { 1, 2, 2, 3 };

    for ( size_t i = 0; i < N; ++i )
    {
        if ( try_basis(V,n,VB,VCs[i],NCs[i]) )
            return;
    }
    throw std::runtime_error { "Error while computing bounding sphere. No basis found. Case n = 2." };
}

void basis_computation_3( std::array<sphere,4> &V, size_t &n, sphere &VB, sphere B )
{
    // n = 3.
    constexpr size_t N {8};
    const std::array<sphere,4> VCs[N]
    {
        { B },
        { V[0], B },
        { V[1], B },
        { V[2], B },
        { V[0], V[1], B },
        { V[0], V[2], B },
        { V[1], V[2], B },
        { V[0], V[1], V[2], B }
    };
    const size_t NCs[N] { 1, 2, 2, 2, 3, 3, 3, 4 };

    for ( size_t i = 0; i < N; ++i )
    {
        if ( try_basis(V,n,VB,VCs[i],NCs[i]) )
            return;
    }
    throw std::runtime_error { "Error while computing bounding sphere. No basis found. Case n = 3." };
}

void basis_computation_4( std::array<sphere,4> &V, size_t &n, sphere &VB, sphere B )
{
    // n = 4.
    constexpr size_t N {15};
    const std::array<sphere,4> VCs[N]
    {
        { B },
        { V[0], B },
        { V[1], B },
        { V[2], B },
        { V[3], B },
        { V[0], V[1], B },
        { V[0], V[2], B },
        { V[0], V[3], B },
        { V[1], V[2], B },
        { V[1], V[3], B },
        { V[2], V[3], B },
        { V[1], V[2], V[3], B },
        { V[0], V[2], V[3], B },
        { V[0], V[1], V[3], B },
        { V[0], V[1], V[2], B }
    };
    const size_t NCs[N] { 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4 };

    for ( size_t i = 0; i < N; ++i )
    {
        if ( try_basis(V,n,VB,VCs[i],NCs[i]) )
            return;
    }
    throw std::runtime_error { "Error while computing bounding sphere. No basis found. Case n = 4." };
}

/*!
 * \brief Computes a new basis from an old one and a violating ball.
 * \param V  The old basis.
 * \param n  Number of balls in the old basis.
 * \param VB The bounding ball corresponding to the old basis.
 * \param B  The ball that violates the basis, i.e., B is not a subset of VB.
 * \see Kaspar Fischer and Bernd Gärtner, "The Smallest Enclosing Ball of Balls:
 *      Combinatorial Structure and Algorithms", section 3.2.
 *
 * Given an old basis V and its corresponding bounding ball VB as well as a 
 * ball B which is not a subset of VB. This method computes a new basis,
 * whose corresponding bounding ball is a superset of all the balls in the
 * old basis and B. It does so using the brute-force approach, as described
 * by Fischer and Gärtner, i.e., it simply tries all possible combinations
 * in increasing order of size. On exit, V, n, and VB will be overwritten
 * with the new basis.
 */
void basis_computation( std::array<sphere,4> &V, size_t &n, sphere &VB, sphere B )
{
    switch ( n )
    {
    case 1: basis_computation_1( V, n, VB, B ); break;
    case 2: basis_computation_2( V, n, VB, B ); break;
    case 3: basis_computation_3( V, n, VB, B ); break;
    case 4: basis_computation_4( V, n, VB, B ); break;
    }
}

}

/*!
 * \brief Computes the minimal bounding sphere of a set of spheres.
 * \see Kaspar Fischer and Bernd Gärtner, "The Smallest Enclosing Ball of Balls:
 *      Combinatorial Structure and Algorithms", section 3.2.
 * \return The minimal bounding sphere of the given set of spheres.
 * \warning Changes the order of the input points.
 */
sphere bounding_sphere( sphere *begin, sphere *end )
{
    if ( begin == end )
        throw std::logic_error { "Cannot compute bounding sphere of an empty set." };

    std::array<sphere,4> V; V[0] = *begin;
    size_t n { 1 }; sphere VB { V[0] };

    while ( ++begin != end )
    {
        // Find furthest sphere.
        sphere *max_sphere { begin };
        real    max_dist   { 0 };
        for ( sphere *i = begin; i != end; ++i )
        {
            real dist = (i->centre - VB.centre).r() + i->radius;
            if ( dist > max_dist )
            {
                max_dist   = dist;
                max_sphere = i;
            }
        }
  
        if ( VB.approx_contains(*max_sphere) ) 
            break;

        basis_computation( V, n, VB, *max_sphere );
        std::swap( *max_sphere, *begin );
    }

    VB.radius *= 1 + sphere::TOL;
    return VB; 
}

/*!
 * \brief Computes the minimal bounding sphere of a set of points.
 * \see Kaspar Fischer and Bernd Gärtner, "The Smallest Enclosing Ball of Balls:
 *      Combinatorial Structure and Algorithms", section 3.2.
 * \return The minimal bounding sphere of the given set of points.
 * \warning Changes the order of the input points.
 */
sphere bounding_sphere( point *begin, point *end )
{
    if ( begin == end )
        throw std::logic_error { "Cannot compute bounding sphere of an empty set." };

    std::array<sphere,4> V;
    V[0].centre = *begin; V[0].radius = 0;
    size_t n { 1 }; sphere VB { V[0] };

    while ( ++begin != end )
    {
        // Find furthest point.
        point *max_point { begin };
        real   max_dist  { 0 };
        for ( point *i = begin; i != end; ++i )
        {
            real dist = (*i - VB.centre).r();
            if ( dist > max_dist )
            {
                max_dist  = dist;
                max_point = i;
            }
        }
  
        if ( VB.approx_contains(*max_point) ) 
            break;

        basis_computation( V, n, VB, sphere { *max_point, 0 } );
        std::swap( *max_point, *begin );
    }

    VB.radius *= 1 + sphere::TOL;
    return VB; 
}

}

