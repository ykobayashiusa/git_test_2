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

namespace bem
{

/*!
 * \brief Construct from given barycentric coordinates.
 */
constexpr bary::bary( real p_z0, real p_z1, real p_z2 ):
 z0 { p_z0 }, z1 { p_z1 }, z2{ p_z2 }
{}

/*!
 * \brief Construct from local cartesian coordinates.
 */
constexpr bary::bary( real xi, real eta ):
 z0 { 1 - xi }, z1 { xi - eta }, z2 { eta }
{}

constexpr bary operator+( const bary& p1, const bary& p2 )
{
	return bary( p1.z0 + p2.z0,
	             p1.z1 + p2.z1,
	             p1.z2 + p2.z2 );
}

constexpr bary operator-( const bary& p1, const bary& p2 )
{
	return bary( p1.z0 - p2.z0,
	             p1.z1 - p2.z1,
	             p1.z2 - p2.z2 );
}

constexpr bary operator*( const real  f, const bary& p )
{
	return bary( f*p.z0, f*p.z1, f*p.z2 );
}

constexpr bary operator*( const bary& p, const real  f )
{
	return bary( f*p.z0, f*p.z1, f*p.z2 );
}

constexpr real bary::xi() const
{
	return z1 + z2;
}

constexpr real bary::eta() const
{
	return z2;
}

}

