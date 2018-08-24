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

namespace geometry
{

constexpr bary2d::bary2d( real p_z0, real p_z1, real p_z2 ) noexcept:
 z0 { p_z0 }, z1 { p_z1 }, z2{ p_z2 }
{}

constexpr bary2d::bary2d( real xi, real eta ) noexcept:
 z0 { 1 - xi }, z1 { xi - eta }, z2 { eta }
{}

constexpr bary2d operator+( const bary2d& p1, const bary2d& p2 ) noexcept
{
	return bary2d( p1.z0 + p2.z0,
	               p1.z1 + p2.z1,
	               p1.z2 + p2.z2 );
}

constexpr bary2d operator-( const bary2d& p1, const bary2d& p2 ) noexcept
{
	return bary2d( p1.z0 - p2.z0,
	               p1.z1 - p2.z1,
	               p1.z2 - p2.z2 );
}

constexpr bary2d operator*( const real  f, const bary2d& p ) noexcept
{
	return bary2d( f*p.z0, f*p.z1, f*p.z2 );
}

constexpr bary2d operator*( const bary2d& p, const real  f ) noexcept
{
	return bary2d( f*p.z0, f*p.z1, f*p.z2 );
}

constexpr bary2d operator/( const bary2d& p, const real  f ) noexcept
{
	return bary2d( p.z0/f, p.z1/f, p.z2/f );
}

constexpr real bary2d::xi() const noexcept
{
	return z1 + z2;
}

constexpr real bary2d::eta() const noexcept
{
	return z2;
}







constexpr bary3d::bary3d( real p_z0, real p_z1, real p_z2, real p_z3 ) noexcept:
 z0 { p_z0 }, z1 { p_z1 }, z2 { p_z2 }, z3 { p_z3 }
{}

constexpr bary3d::bary3d( real xi, real eta, real zeta ) noexcept:
 z0 { 1 - xi - eta - zeta }, z1 { xi }, z2 { eta }, z3 { zeta }
{}

constexpr bary3d operator+( const bary3d& p1, const bary3d& p2 ) noexcept
{
	return bary3d( p1.z0 + p2.z0,
	               p1.z1 + p2.z1,
	               p1.z2 + p2.z2,
                   p1.z3 + p2.z3 );
}

constexpr bary3d operator-( const bary3d& p1, const bary3d& p2 ) noexcept
{
	return bary3d( p1.z0 - p2.z0,
	               p1.z1 - p2.z1,
	               p1.z2 - p2.z2,
                   p1.z3 - p2.z3 );
}

constexpr bary3d operator*( const real  f, const bary3d& p ) noexcept
{
	return bary3d( f*p.z0, f*p.z1, f*p.z2, f*p.z3 );
}

constexpr bary3d operator*( const bary3d& p, const real  f ) noexcept
{
	return bary3d( f*p.z0, f*p.z1, f*p.z2, f*p.z3 );
}

constexpr bary3d operator/( const bary3d& p, const real  f ) noexcept
{
	return bary3d( p.z0/f, p.z1/f, p.z2/f, p.z3/f );
}

constexpr real bary3d::xi() const noexcept
{
	return z1;
}

constexpr real bary3d::eta() const noexcept
{
	return z2;
}

constexpr real bary3d::zeta() const noexcept
{
    return z3;
}

}

namespace std
{

template <>
struct hash<geometry::bary2d>
{
    size_t operator()( geometry::bary2d pos ) const noexcept
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


template <>
struct hash<geometry::bary3d>
{
    size_t operator()( geometry::bary3d pos ) const noexcept
    {
        using boost::hash_combine;
        std::hash<decltype(pos.z0)> coord_hash;

        size_t result { 0 };
        hash_combine( result, coord_hash( pos.z0 ) );
        hash_combine( result, coord_hash( pos.z1 ) );
        hash_combine( result, coord_hash( pos.z2 ) );
        hash_combine( result, coord_hash( pos.z3 ) );

        return result;
    } 
};

}

