/*
 * Copyright (C) 2014 Matthias Kirchhart
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
#ifndef BEM_GEOID_H
#define BEM_GEOID_H

#include <limits>
#include <functional>

#include <boost/functional/hash.hpp> // For hash_combine.

#include "geometry/point.h"

namespace bem
{

struct geoid
{
    unsigned char   dim;
    unsigned char   lvl; 
    geometry::point pos;

    // The default constructor now creates NO_GEOID by default.
    constexpr geoid() noexcept:
    dim { std::numeric_limits<unsigned char>::max() },
    lvl { std::numeric_limits<unsigned char>::max() },
    pos { std::numeric_limits<real>::max(),
          std::numeric_limits<real>::max(),
          std::numeric_limits<real>::max() }
    {}

    // This constructor allows us to emulate aggregate initialisation.
    constexpr geoid( unsigned char p_dim, unsigned char p_lvl, geometry::point p_pos ):
    dim { p_dim }, lvl { p_lvl }, pos { p_pos }
    {}

};

constexpr bool operator==( const geoid& lhs, const geoid& rhs ) noexcept
{
    return lhs.dim == rhs.dim &&
           lhs.lvl == rhs.lvl &&
           lhs.pos == rhs.pos;
}

constexpr bool operator!=( const geoid& lhs, const geoid& rhs ) noexcept
{
    return ! ( lhs == rhs );
}

// A global constant representing no object. Analogous to nullptr.
constexpr geoid NO_GEOID;

}

namespace std
{

template <>
struct hash<bem::geoid>
{
    size_t operator()( bem::geoid id ) const noexcept
    {
        using boost::hash_combine;
        std::hash<unsigned char>      uchar_hash;
        std::hash<decltype(id.pos.x)> coord_hash;

        size_t result { 0 };
        hash_combine( result, uchar_hash( id.dim   ) );
        hash_combine( result, uchar_hash( id.lvl   ) );
        hash_combine( result, coord_hash( id.pos.x ) );
        hash_combine( result, coord_hash( id.pos.y ) );
        hash_combine( result, coord_hash( id.pos.z ) );

        return result;
    } 
};

}

#endif

