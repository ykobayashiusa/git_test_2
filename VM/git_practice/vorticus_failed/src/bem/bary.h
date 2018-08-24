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
#ifndef BEM_BARY_H
#define BEM_BARY_H

#include "types.h"

#include <cmath>
#include <functional>

#include <boost/functional/hash.hpp>

namespace bem
{

/*!
 * \brief Points in the reference triangle using barycentric coordinates.
 */
struct bary
{
	constexpr bary() = default;
	constexpr bary( const bary& ) = default;
	bary& operator=( const bary& ) = default;

	constexpr bary( real p_z0, real p_z1, real p_z2 );
	constexpr bary( real xi, real eta );

	constexpr real  xi() const;
	constexpr real eta() const;

	real z0 {0};
	real z1 {0};
	real z2 {0};
};

constexpr bary operator+( const bary& p1, const bary& p2 );
constexpr bary operator-( const bary& p1, const bary& p2 );

constexpr bary operator*( const real  f, const bary& p );
constexpr bary operator*( const bary& p, const real  f );

}

namespace std
{

template <>
struct hash<bem::bary>
{
    size_t operator()( bem::bary pos ) const noexcept
    {
        using boost::hash_combine;
        std::hash<decltype(pos.z0)> coord_hash;

        size_t result { 0 };
        hash_combine( result, coord_hash( pos.z0 ) );
        hash_combine( result, coord_hash( pos.z1 ) );
        hash_combine( result, coord_hash( pos.z2 ) );

        return result;
    } 
};

}

#include "bary.tpp"
#endif

