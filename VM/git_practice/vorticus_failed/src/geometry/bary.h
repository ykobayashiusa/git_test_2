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
#ifndef GEOMETRY_BARY_H
#define GEOMETRY_BARY_H

#include "types.h"

#include <cmath>
#include <functional>

#include <boost/functional/hash.hpp>

namespace geometry
{

struct bary2d
{
	constexpr bary2d() noexcept = default;
	constexpr bary2d( const bary2d& ) noexcept = default;
	bary2d& operator=( const bary2d& ) noexcept = default;

	constexpr bary2d( real p_z0, real p_z1, real p_z2 ) noexcept;
	constexpr bary2d( real xi, real eta ) noexcept;

	constexpr real   xi() const noexcept;
	constexpr real  eta() const noexcept;

	real z0 {0};
	real z1 {0};
	real z2 {0};
};

constexpr bary2d operator+( const bary2d& p1, const bary2d& p2 ) noexcept;
constexpr bary2d operator-( const bary2d& p1, const bary2d& p2 ) noexcept;

constexpr bary2d operator*( const real  f, const bary2d& p ) noexcept;
constexpr bary2d operator*( const bary2d& p, const real  f ) noexcept;
constexpr bary2d operator/( const bary2d& p, const real  f ) noexcept;



struct bary3d
{
	constexpr bary3d() noexcept = default;
	constexpr bary3d( const bary3d& ) noexcept = default;
	bary3d& operator=( const bary3d& ) noexcept = default;

	constexpr bary3d( real p_z0, real p_z1, real p_z2, real p_z3 ) noexcept;
	constexpr bary3d( real xi, real eta, real zeta ) noexcept;

	constexpr real   xi() const noexcept;
	constexpr real  eta() const noexcept;
	constexpr real zeta() const noexcept;

	real z0 {0};
	real z1 {0};
	real z2 {0};
	real z3 {0};
};

constexpr bary3d operator+( const bary3d& p1, const bary3d& p2 ) noexcept;
constexpr bary3d operator-( const bary3d& p1, const bary3d& p2 ) noexcept;

constexpr bary3d operator*( const real    f, const bary3d& p ) noexcept;
constexpr bary3d operator*( const bary3d& p, const real    f ) noexcept;
constexpr bary3d operator/( const bary3d& p, const real    f ) noexcept;

}

#include "geometry/bary.tpp"
#endif

