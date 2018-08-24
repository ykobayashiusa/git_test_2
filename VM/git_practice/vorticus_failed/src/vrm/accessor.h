/*
 * Copyright (C) 2015 Matthias Kirchhart
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
#ifndef VRM_ACCESSOR_H
#define VRM_ACCESSOR_H

#include "vrm/particle.h"

namespace vrm
{

template <uint stages>
class accessor;

template <>
class accessor<1>
{
public:
    using point = geometry::point;

    point operator()( const particle<1> &p ) const noexcept
    {
        return p.x;
    }

    point x( const particle<1> &p ) const noexcept
    {
        return p.x;
    }

    point G( const particle<1> &p ) const noexcept
    {
        return p.G;
    }

    point  u( const particle<1> &p ) const noexcept
    {
        return p.k_x[ 0 ];
    }

    point& u( particle<1> &p )
    {
        return p.k_x[ 0 ];
    }

    point  dG( const particle<1> &p ) const noexcept
    {
        return p.k_G[ 0 ];
    }

    point& dG( particle<1> &p )
    {
        return p.k_G[ 0 ];
    }

    real dt() const noexcept
    {
        return m_dt;
    }

    void dt( real val ) noexcept
    {
        assert( val > 0 );
        m_dt = val;
    }

    uint stage() const noexcept
    {
        return 0;
    }

    void stage( uint )
    {
    }

    void advance( particle<1> &p )
    {
        p.x += m_dt * p.k_x[0];
        p.G += m_dt * p.k_G[0];
        p.k_x[0] = point {};
        p.k_G[0] = point {};
    }

private:
    real m_dt     { 1 };
};

template <>
class accessor<2>
{
public:
    using point = geometry::point;

    point operator()( const particle<2> &p ) const noexcept
    {
        switch ( stage() )
        {
        case 0: return p.x;
        case 1: return p.x + m_dt*p.k_x[ 0 ];
        default: return p.x;
        }
    }

    point x( const particle<2> &p ) const noexcept
    {
        switch ( stage() )
        {
        case 0: return p.x;
        case 1: return p.x + m_dt*p.k_x[ 0 ];
        default: return p.x;
        }
    }

    point G( const particle<2> &p ) const noexcept
    {
        switch ( stage() )
        {
        case 0: return p.G;
        case 1: return p.G + m_dt*p.k_G[ 0 ];
        default: return p.G;
        }
    }

    point  u( const particle<2> &p ) const noexcept
    {
        return p.k_x[ stage() ];
    }

    point& u( particle<2> &p )
    {
        return p.k_x[ stage() ];
    }

    point dG( const particle<2> &p ) const noexcept
    {
        return p.k_G[ stage() ];
    }

    point& dG( particle<2> &p )
    {
        return p.k_G[ stage() ];
    }

    real dt() const noexcept
    {
        return m_dt;
    }

    void dt( real val ) noexcept
    {
        assert( val > 0 );
        m_dt = val;
    }

    uint stage() const noexcept
    {
        return m_stage;
    }

    void stage( uint s )
    {
        assert( s <= 1 );
        m_stage = s;
    }

    void advance( particle<2> &p )
    {
        p.x += m_dt * ( .5*p.k_x[0] + .5*p.k_x[1] );
        p.G += m_dt * ( .5*p.k_G[0] + .5*p.k_G[1] );
        p.k_x[0] = p.k_x[1] = point {};
        p.k_G[0] = p.k_G[1] = point {};
        stage(0);
    }

private:
    real m_dt     { 1 };
    uint m_stage  { 0 };
};

template <>
class accessor<4>
{
public:
    using point = geometry::point;

    point operator()( const particle<4> &p ) const noexcept
    {
        switch ( stage() )
        {
        case 0: return p.x;
        case 1: return p.x + 0.5*m_dt*p.k_x[ 0 ];
        case 2: return p.x + 0.5*m_dt*p.k_x[ 1 ];
        case 3: return p.x +     m_dt*p.k_x[ 2 ];
        default: return p.x;
        }
    }

    point x( const particle<4> &p ) const noexcept
    {
        switch ( stage() )
        {
        case 0: return p.x;
        case 1: return p.x + 0.5*m_dt*p.k_x[ 0 ];
        case 2: return p.x + 0.5*m_dt*p.k_x[ 1 ];
        case 3: return p.x +     m_dt*p.k_x[ 2 ];
        default: return p.x;
        }
    }

    point G( const particle<4> &p ) const noexcept
    {
        switch ( stage() )
        {
        case 0: return p.G;
        case 1: return p.G + 0.5*m_dt*p.k_G[ 0 ];
        case 2: return p.G + 0.5*m_dt*p.k_G[ 1 ];
        case 3: return p.G +     m_dt*p.k_G[ 2 ];
        default: return p.G;
        }
    }

    point  u( const particle<4> &p ) const noexcept
    {
        return p.k_x[ stage() ];
    }

    point& u( particle<4> &p )
    {
        return p.k_x[ stage() ];
    }

    point  dG( const particle<4> &p ) const noexcept
    {
        return p.k_G[ stage() ];
    }

    point& dG( particle<4> &p )
    {
        return p.k_G[ stage() ];
    }

    real dt() const noexcept
    {
        return m_dt;
    }

    void dt( real val ) noexcept
    {
        assert( val > 0 );
        m_dt = val;
    }

    uint stage() const noexcept
    {
        return m_stage;
    }

    void stage( uint s )
    {
        assert( s <= 3 );
        m_stage = s;
    }

    void advance( particle<4> &p )
    {
        p.x += m_dt * ( p.k_x[0]/6. + p.k_x[1]/3. + p.k_x[2]/3. + p.k_x[3]/6. );
        p.G += m_dt * ( p.k_G[0]/6. + p.k_G[1]/3. + p.k_G[2]/3. + p.k_G[3]/6. );
        p.k_x[0] = p.k_x[1] = p.k_x[2] = p.k_x[3] = point {};
        p.k_G[0] = p.k_G[1] = p.k_G[2] = p.k_G[3] = point {};
        stage(0);
    }

private:
    real m_dt     { 1 };
    uint m_stage  { 0 };
};

}

#endif

