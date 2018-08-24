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

inline
box bounding_box( const box &p, const box &q ) noexcept
{
    using std::min;
    using std::max;

    return box
    {
        point { min(p.min.x,q.min.x), min(p.min.y, q.min.y), min(p.min.z,q.min.z) },
        point { max(p.max.x,q.max.x), max(p.max.y, q.max.y), max(p.max.z,q.max.z) }
    };
}

inline
box bounding_box( const box &p, const point &q ) noexcept
{
    using std::min;
    using std::max;

    return box
    {
        point { min(p.min.x, q.x), min(p.min.y, q.y), min(p.min.z, q.z) },
        point { max(p.max.x, q.x), max(p.max.y, q.y), max(p.max.z, q.z) }
    };
}

inline
box bounding_box( const point &p, const box &q ) noexcept
{
    return bounding_box(q,p);
}

inline
real dist( const box &p, const box &q ) noexcept
{
    using std::min;
    using std::max;
    point dist {};
    dist.x = p.centre().x < q.centre().x ? max( 0., q.min.x - p.max.x ) : max( 0., p.min.x - q.max.x );
    dist.y = p.centre().y < q.centre().y ? max( 0., q.min.y - p.max.y ) : max( 0., p.min.y - q.max.y );
    dist.z = p.centre().z < q.centre().z ? max( 0., q.min.z - p.max.z ) : max( 0., p.min.z - q.max.z );
    return dist.r();
}


// Returns the distance of the point in b which is closest to p.
inline
real min_dist( const box &b, const point &p ) noexcept
{
    using std::abs;

    point q;
    if ( b.min.x <= p.x && p.x <= b.max.x ) q.x = p.x;
    else q.x = ( abs(b.min.x - p.x) < abs(b.max.x-p.x) ) ? b.min.x : b.max.x;

    if ( b.min.y <= p.y && p.y <= b.max.x ) q.y = p.y;
    else q.y = ( abs(b.min.y - p.y) < abs(b.max.y-p.y) ) ? b.min.y : b.max.y;

    if ( b.min.z <= p.z && p.z <= b.max.z ) q.z = p.z;
    else q.z = ( abs(b.min.z - p.z) < abs(b.max.z-p.z) ) ? b.min.z : b.max.z;

    return (p - q).r();
}


// Returns the distance of the point in b which is furthest from p.
inline
real max_dist( const box &b, const point &p ) noexcept
{
    using std::abs;
    point q;

    q.x = ( abs(b.min.x - p.x) > abs(b.max.x - p.x) ) ? b.min.x : b.max.x;
    q.y = ( abs(b.min.y - p.y) > abs(b.max.y - p.y) ) ? b.min.y : b.max.y;
    q.z = ( abs(b.min.z - p.z) > abs(b.max.z - p.z) ) ? b.min.z : b.max.z;

    return (p - q).r();
}

}

