/*
 * Copyright (C) 2016 Matthias Kirchhart
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
#include <parallel/algorithm>

namespace vrm
{

namespace laplacian
{

constexpr real  cone_angle { math::pi/4 };
constexpr geometry::point cone_dir[8]
{
    geometry::point {  1                 ,  0                 , 0 },
    geometry::point {  0.7071067811865475,  0.7071067811865475, 0 },
    geometry::point {  0                 ,  1                 , 0 },
    geometry::point { -0.7071067811865475,  0.7071067811865475, 0 },
    geometry::point { -1                 ,  0                 , 0 },
    geometry::point { -0.7071067811865475, -0.7071067811865475, 0 },
    geometry::point {  0                 , -1                 , 0 },
    geometry::point {  0.7071067811865475, -0.7071067811865475, 0 }
};


template <typename iterator, typename accessor>
auto active_set( iterator begin, iterator end, accessor &get,
                 real h, real C_diff ) -> std::vector<iterator>
{
    std::vector<iterator> vec( end - begin );
    real sum { 0 };
    #pragma omp parallel for reduction(+:sum)
    for ( decltype(vec.size()) i = 0; i < vec.size(); ++i )
    {
        vec[ i ] = begin + i;
        sum += length(get.G(*vec[i]));
    }

    // Sort descending according to circulation.
    // std::sort( vec.begin(), vec.end(),
    //            [&get]( iterator p1, iterator p2 )
    //            { return length(get.G(*p1)) > length(get.G(*p2)); } );
    __gnu_parallel::sort( vec.begin(), vec.end(),
                          [&get]( iterator p1, iterator p2 )
                          { return length(get.G(*p1)) > length(get.G(*p2)); } );

    const real h3 { h*h*h };
    real excluded { 0 };
    while ( vec.size() )
    {
        excluded += length(get.G(*vec.back()));
        if ( excluded <= C_diff*h3*sum ) vec.pop_back();
        else break;
    }

    // To preserve spatial locality: sort according to original order, which
    // should be close to a Morton z-order curve.
    //std::sort( vec.begin(), vec.end() );
    __gnu_parallel::sort( vec.begin(), vec.end() );


    return std::move(vec);
}

template <typename iterator, typename accessor>
bool
small_neighbourhood_2d( iterator it, const std::vector<iterator> &full,
                        accessor &get, std::array<iterator,8> &result )
{
    using geometry::point;
    using size_type = typename std::vector<iterator>::size_type;

    const iterator nil {};
    result.fill(nil);

    const point x_p = get(*it);
    size_t found_count { 0 };


    // For each particle: find the corresponding segment it belongs to.
    // This can be done by checking the angle between the particle and
    // the segment's centre line. It should be less then half the segments
    // opening angle, so less then 45/2 = 22.5 degrees. We make use of
    // arccos' monotonicity and avoid using the trigonometric functions.
    constexpr real cos225 = .9238795325112867561;
    for ( size_type k = 0; k < full.size() && found_count < 8; ++k )
    {
        point r_kp = get(*full[k]) - x_p;
        real  r    = r_kp.r();

        for ( size_t j = 0; j < 8; ++j )
        {
            //real cos_phi = scal_prod( r_kp, cone_dir[ j ] ) / r;

            geometry::point conedir{ std::cos( (math::pi/180.)*10 + j*(math::pi/4.) ),
                                     std::sin( (math::pi/180.)*10 + j*(math::pi/4.) ), 0 };
            real cos_phi = scal_prod( r_kp, conedir  ) / r;
            if ( cos_phi >= cos225 )
            {
                if ( result[ j ] == nil )
                {
                    ++found_count;
                    result[ j ] = full[ k ];
                }
                else if ( r*r < ( get(*result[j]) - x_p ).r2() )
                {
                    // Do not increase found_count!
                    result[ j ] = full[ k ];
                }
                break; // No need to check the other segments.
            }
        }
    }

    return found_count == 8;
}

template <typename iterator, typename accessor, typename finder>
void full_neighbourhood( iterator p, accessor &get, finder &find,
                         real h, std::vector<iterator> &v )
{
    v.clear();
    const geometry::point centre = get(*p);
    auto pred = [&get,centre,h]( iterator it ) -> bool
    {
        return ( centre - get(*it) ).r2() <    4*h*h &&
               ( centre - get(*it) ).r2() > 0.25*h*h;
    };
    geometry::box b( geometry::point { centre.x - 2*h, centre.y - 2*h, centre.z - 2*h },
                     geometry::point { centre.x + 2*h, centre.y + 2*h, centre.z + 2*h } );
    find.query( b, pred, std::back_inserter( v ) );
}

                         
template <typename iterator, typename accessor,
          typename finder,   typename maker>
std::vector<std::array<iterator,8>>
create_neighbourhoods_2d( const std::vector<iterator> &particles,
                          accessor &get, finder &find, maker &make,
                          real h )
{
    using size_type = typename std::vector<iterator>::size_type;
    std::vector<std::array<iterator,8>> result( particles.size() );

    std::mutex m {};
    std::vector<size_type> failed;
    std::vector<iterator>  full;

    iterator nil {};

    // This loop does the main work in the process of finding the particle
    // neighbourhoods. The full-vector is firstprivate, such that every thread
    // has its own copy of it, which acts as a buffer, avoiding reallocation
    // in every iteration.
    #pragma omp parallel for firstprivate(full) schedule(dynamic)
    for ( size_type i = 0; i < particles.size(); ++i )
    {
        iterator p = particles[ i ];
        full_neighbourhood( p, get, find, h, full );
        bool success = small_neighbourhood_2d( p, full, get, result[ i ] );
        if ( success == false )
        {
            std::unique_lock<std::mutex> l(m);
            failed.push_back( i );
        }
    }

    // This loop fills holes in the particles’ neighbourhoods. Due to the
    // write access, it cannot be (easily?) parallelised.
    for ( size_type i: failed )
    {
        iterator p = particles[ i ];
        full_neighbourhood( p, get, find, h, full );
        small_neighbourhood_2d( p, full, get, result[ i ] );
        for ( size_type j = 0; j < 8; ++j )
        {
            if ( result[ i ][ j ] == nil )
            {
                result[ i ][ j ] = make( get(*p) + 1.5*h*cone_dir[j] );
            }
        }
    }

    return std::move(result);
}

template <typename iterator, typename accessor>
bool small_solution_2d( iterator p, const std::array<iterator,8> &neigh,
                        accessor &get, real h,
                        typename stiffness_matrix<iterator>::iterator target )
{
    arma::vec::fixed<5>   rhs { 0, 0, 2., 2., 0 };
    arma::vec::fixed<8>   f_pq;
    arma::mat::fixed<5,8> V;

    for ( size_t q_idx = 0; q_idx < 8; ++q_idx )
    {
        iterator q = neigh[ q_idx ];
        const geometry::point x_pq = (get(*q) - get(*p))/h;
        V( 0, q_idx ) = x_pq.x;
        V( 1, q_idx ) = x_pq.y;
        
        V( 2, q_idx ) = x_pq.x * x_pq.x;
        V( 3, q_idx ) = x_pq.y * x_pq.y;
        V( 4, q_idx ) = x_pq.x * x_pq.y;
    }

    int ret = math::rev_nnsolve( 5, 8, V.memptr(), rhs.memptr(), f_pq.memptr() );
    //int ret = math::nnsolve( V, rhs, f_pq );
    if ( ret != 0 )
    {
        std::cerr << "Instability in Simplex solver.\n";
        return false;
    }
    if ( norm( V*f_pq - rhs, 2 ) > 1e-9 )
    {
        // No feasible solution.
        std::cout << "No feasible solution! Residual: " << norm( V*f_pq - rhs, 2 ) << ".\n";
        return false;
    }

    // Rescale.
    f_pq /= h*h;

    auto end = target + 6;
    std::get<1>(*target) =  p;
    std::get<2>(*target) = -sum(f_pq);
    ++target;
    for ( size_t q_idx = 0; q_idx < 8 && target != end; ++q_idx )
    {
        if ( f_pq(q_idx) != 0 )
        {
            std::get<1>(*target) = neigh[q_idx];
            std::get<2>(*target) = f_pq(q_idx);
            ++target;
        }
    }

    return true;
}

template <typename iterator, typename accessor>
bool full_solution_2d( iterator p, const std::vector<iterator> &neigh,
                       accessor &get, real h,
                       typename stiffness_matrix<iterator>::iterator target )
{
    math::vector   rhs { 0, 0, 2., 2., 0 };
    math::vector   f_pq( neigh.size() );
    math::matrix   V(5,neigh.size());

    using size_type = decltype(neigh.size());
    for ( size_type q_idx = 0; q_idx < neigh.size(); ++q_idx )
    {
        iterator q = neigh[ q_idx ];
        const geometry::point x_pq = (get(*q) - get(*p))/h;
        V( 0, q_idx ) = x_pq.x;
        V( 1, q_idx ) = x_pq.y;
        
        V( 2, q_idx ) = x_pq.x * x_pq.x;
        V( 3, q_idx ) = x_pq.y * x_pq.y;
        V( 4, q_idx ) = x_pq.x * x_pq.y;
    }

    //int ret = math::nnsolve( V, rhs, f_pq );
    int ret = math::rev_nnsolve( 5, neigh.size(), V.memptr(), rhs.memptr(), f_pq.memptr() );
    if ( ret != 0 )
    {
        std::cerr << "Instability in Simplex solver.\n";
        return false;
    }
    if ( norm( V*f_pq - rhs, 2 ) > 1e-9 )
    {
        // No feasible solution.
        std::cout << "No feasible solution! Residual: " << norm( V*f_pq - rhs, 2 ) << ".\n";
        return false;
    }

    // Rescale
    f_pq /= h*h;

    auto end = target + 6;
    std::get<1>(*target) =  p;
    std::get<2>(*target) = -sum(f_pq);
    ++target;
    for ( size_t q_idx = 0; q_idx < neigh.size() && target != end; ++q_idx )
    {
        if ( f_pq(q_idx) != 0 )
        {
            std::get<1>(*target) = neigh[q_idx];
            std::get<2>(*target) = f_pq(q_idx);
            ++target;
        }
    }

    return true;
}

template <typename iterator, typename accessor,
          typename finder,   typename maker>
auto F_2d( iterator begin, iterator end,
           accessor &get, finder &find, maker &make,
           real h, real C_diff )
-> stiffness_matrix<iterator>
{
    const std::vector<iterator>               diffused   = active_set( begin, end, get, h, C_diff );
    std::vector<std::array<iterator,8>>       neighbours = create_neighbourhoods_2d( diffused, get, find, make, h );

    stiffness_matrix<iterator> result( diffused.size() * 6 );
    using size_type = decltype( diffused.size() );
    #pragma omp parallel for
    for ( size_type i = 0; i < diffused.size(); ++i )
    {
        std::get<0>(result[ 6*i + 0 ]) = diffused[ i ];
        std::get<0>(result[ 6*i + 1 ]) = diffused[ i ];
        std::get<0>(result[ 6*i + 2 ]) = diffused[ i ];
        std::get<0>(result[ 6*i + 3 ]) = diffused[ i ];
        std::get<0>(result[ 6*i + 4 ]) = diffused[ i ];
        std::get<0>(result[ 6*i + 5 ]) = diffused[ i ];
    }

    std::vector<iterator> full;
    #pragma omp parallel for schedule(dynamic) firstprivate(full)
    for ( size_type i = 0; i < diffused.size(); ++i )
    {
        bool success = small_solution_2d( diffused[ i ], neighbours[ i ], get, h, result.begin() + 6*i );
        if ( success == false )
        {
            full_neighbourhood( diffused[ i ], get, find, h, full );
            success = full_solution_2d( diffused[ i ], full, get, h, result.begin() + 6*i );
            if ( success == false )
            {
                std::cerr << "No feasible solution for moment equations!\n";
                std::exit(-1);
            }
        }
    }

    return std::move(result);
}

}


template <typename iterator, typename accessor,
          typename finder,   typename maker>
auto laplacian_2d( iterator begin, iterator end,
                   accessor &get, finder &find, maker &make,
                   real h, real C_diff )
-> stiffness_matrix<iterator>
{
    return laplacian::F_2d( begin, end, get, find, make, h, C_diff );
}

}

