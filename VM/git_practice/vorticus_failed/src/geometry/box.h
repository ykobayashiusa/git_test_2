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
#ifndef GEOMETRY_BOX_H
#define GEOMETRY_BOX_H

#include "geometry/point.h"

namespace geometry
{

struct box
{
    constexpr box() noexcept = default;
    constexpr box( point minp, point maxp ) noexcept:
    min { minp }, max { maxp } {}

    constexpr point centre() const noexcept
    { return (min+max)/2; }

    real diameter() const noexcept
    { return (max-min).r(); }
    
    constexpr real diameter2() const noexcept
    { return (max-min).r2(); }

    constexpr point diag() const noexcept
    { return max - min; }

    point min {};
    point max {};
};

box bounding_box( const box   &p, const box   &q ) noexcept;
box bounding_box( const box   &p, const point &q ) noexcept;
box bounding_box( const point &p, const box   &q ) noexcept;

real dist( const box &p, const box &q ) noexcept;
real min_dist( const box &b, const point &p ) noexcept;
real max_dist( const box &b, const point &p ) noexcept;

bool intersects_triangle   ( const box &B, const std::array<point,3> &T ) noexcept;
bool intersects_tetrahedron( const box &B, const std::array<point,4> &T ) noexcept;

}

#include "geometry/box.tpp"
#endif

