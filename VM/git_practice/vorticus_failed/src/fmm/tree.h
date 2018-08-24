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
#ifndef FMM_TREE_H
#define FMM_TREE_H

#include "geometry/point.h"
#include "geometry/sphere.h"

#include <cmath>
#include <memory>
#include <algorithm>

namespace fmm
{

template <typename coeff_t, typename iterator>
struct tree
{
    tree() = delete;
    tree( const tree&  ) = delete;
    tree(       tree&& ) = default;
    tree& operator=( const tree&  ) = delete;
    tree& operator=(       tree&& ) = default;
   ~tree() = default;

    template <typename strategy>
    tree( iterator begin, iterator end,
          strategy &S, size_t leaf_max );

    struct box;
    std::vector<iterator> elems;
    std::unique_ptr<box>  root;
};

template <typename coeff_t, typename iterator>
struct tree<coeff_t,iterator>::box
{
    using point = geometry::point;

    box() = delete;
    box( const box&  ) = delete;
    box(       box&& ) = delete;
    box& operator=( const box&  ) = delete;
    box& operator=(       box&& ) = delete;


    template <typename strategy> box( iterator *first, iterator *last,
                                      strategy &S, size_t leaf_max );

    template <typename strategy> void find_split_centre( strategy &S );
    template <typename strategy> void resize( strategy &S );
    template <typename strategy> void split( strategy &S, size_t leaf_max );

    bool is_leaf() const;

    real  length {};
    point centre {};

    iterator *begin, *end;
    coeff_t   coeffs {};

    std::unique_ptr<box> children[ 8 ];
};



template <typename coeff_t, typename iterator>
template <typename strategy>
tree<coeff_t,iterator>::tree( iterator begin, iterator end,
                              strategy &S, size_t leaf_max ):
elems (std::distance(begin,end)), root { nullptr }
{
    for ( size_t i = 0; i < elems.size(); ++i )
        elems[ i ] = begin++;

    root.reset( new box( elems.data(), elems.data() + elems.size(), S, leaf_max ) );
}

template <typename coeff_t, typename iterator>
template <typename strategy>
void tree<coeff_t,iterator>::box::split( strategy &S, size_t leaf_max )
{
    if ( end - begin <= (int) leaf_max ||
         end - begin == 1 ) return;

    auto xless = [&S,this]( iterator p ){ return S.pos(p).x < centre.x; };
    auto yless = [&S,this]( iterator p ){ return S.pos(p).y < centre.y; };
    auto zless = [&S,this]( iterator p ){ return S.pos(p).z < centre.z; };

    iterator* oct[9];
    oct[0] = begin;
    oct[8] = end;

    using std::partition;
    oct[4] = partition( oct[0], oct[8], xless );

    oct[2] = partition( oct[0], oct[4], yless );
    oct[1] = partition( oct[0], oct[2], zless );
    oct[3] = partition( oct[2], oct[4], zless );

    oct[6] = partition( oct[4], oct[8], yless );
    oct[5] = partition( oct[4], oct[6], zless );
    oct[7] = partition( oct[6], oct[8], zless );

    for ( uint i = 0; i < 8; ++i )
    {
        if ( oct[i] != oct[i+1] )
        {
            #pragma omp task shared(S)
            children[ i ].reset( new box( oct[i], oct[i+1], S, leaf_max ) );
        }
    }
    #pragma omp taskwait
}


template <typename coeff_t, typename iterator>
template <typename strategy>
void tree<coeff_t,iterator>::box::find_split_centre( strategy &S )
{
    using std::min;
    using std::max;
    using geometry::point;

    point P { S.pos(*begin) };
    point minp { P };
    point maxp { P };
    for ( iterator* p = begin; p != end; ++p )
    {
        P = S.pos(*p);

        minp.x = min( minp.x, P.x );
        minp.y = min( minp.y, P.y );
        minp.z = min( minp.z, P.z );

        maxp.x = max( maxp.x, P.x );
        maxp.y = max( maxp.y, P.y );
        maxp.z = max( maxp.z, P.z );
    }
    centre = (minp+maxp)/2;
}

template <typename coeff_t, typename iterator>
template <typename strategy>
void tree<coeff_t,iterator>::box::resize( strategy &S )
{
    using std::min;
    using std::max;
    using geometry::sphere;

    if ( is_leaf() )
    {
        point minp { S.pos(*begin) };
        point maxp { S.pos(*begin) };
        point min_corner, max_corner;
        for ( iterator* p = begin; p != end; ++p )
        {
            S.bounding_box( *p, min_corner, max_corner );

            minp.x = min( minp.x, min_corner.x );
            minp.y = min( minp.y, min_corner.y );
            minp.z = min( minp.z, min_corner.z );

            maxp.x = max( maxp.x, max_corner.x );
            maxp.y = max( maxp.y, max_corner.y );
            maxp.z = max( maxp.z, max_corner.z );
        }
        centre = (minp + maxp)     / 2;
        length = (maxp - minp).r() / 2;
    }
    else
    {
        std::array<sphere,8> spheres; size_t count { 0 };
        for ( size_t i = 0; i < 8; ++i )
        {
            if ( children[i] )
                spheres[ count++ ] = sphere { children[i]->centre, children[i]->length };
        }
        sphere result = bounding_sphere(spheres.data(),spheres.data()+count);
        centre = result.centre;
        length = result.radius;
    }
}

template <typename coeff_t, typename iterator>
template <typename strategy>
tree<coeff_t,iterator>::box::box( iterator *first, iterator *last,
                                  strategy &S, size_t leaf_max ):
begin {first}, end{last}
{
    find_split_centre( S );
    split( S, leaf_max );
    resize( S );
}

template <typename coeff_t, typename iterator>
bool tree<coeff_t,iterator>::box::is_leaf() const
{
    for ( const auto& v: children )
    {
        if ( v ) return false;
    }
    return true;
}

}

#endif

