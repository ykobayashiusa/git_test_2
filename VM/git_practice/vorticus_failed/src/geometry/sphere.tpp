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

namespace geometry
{

constexpr
bool sphere::contains( point p ) const noexcept
{
    return (p-centre).r2() <= radius*radius;
}

constexpr
bool sphere::approx_contains( point p ) const noexcept
{
    return ((p-centre).r2() - (radius*radius)) <= TOL * (radius*radius);
}

constexpr
bool sphere::contains( sphere b ) const noexcept
{
    return ( (radius - b.radius) < 0 ) ?
    ( false ) 
    :
    ( (centre-b.centre).r2() <= (radius - b.radius)*(radius - b.radius) );
}

constexpr
bool sphere::approx_contains( sphere b ) const noexcept
{
    return ( (radius - b.radius) < 0 ) ?
    ( false )
    :
    ( ((centre-b.centre).r2() - (radius-b.radius)*(radius-b.radius)) <= TOL * (radius*radius) );
}

inline
real distance( sphere A, sphere B ) noexcept
{
    return std::max( 0., (A.centre-B.centre).r() - (A.radius+B.radius) );
}

inline
real distance( sphere A, point p ) noexcept
{
    return std::max( 0., (A.centre-p).r() - A.radius );
}

inline
real distance( point p, sphere A ) noexcept
{
    return std::max( 0., (A.centre-p).r() - A.radius );
}

}

