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

namespace vrm
{

template <uint eorder, typename accessor>
template <typename source_iterator>
void particle_strategy<eorder,accessor>::p2m
(
    Mcoeff_t &M, point xc, source_iterator *const begin,
                           source_iterator *const end
) const noexcept
{
    for ( const source_iterator *i = begin; i != end; ++i )
    {
        source_iterator p = *i;
        std::array<cmplx,eorder> transfer = I( get(*p) - xc );

        for ( uint k = 0; k < eorder; ++k )
        {
            M[k] += transfer[k]*(get.G(*p)).z;
        }
    }
}

template <uint eorder, typename accessor>
template <typename source_iterator>
void particle_strategy<eorder,accessor>::p2l
(
    Lcoeff_t &L, point xc, source_iterator *const begin,
                           source_iterator *const end
) const noexcept
{
    for ( const source_iterator *i = begin; i != end; ++i )
    {
        source_iterator p = *i;
        std::array<cmplx,eorder> transfer = O( get(*p) - xc );

        for ( uint k = 0; k < eorder; ++k )
        {
            L[k] += transfer[k]*(get.G(*p)).z;
        }
    }
}

template <uint eorder, typename accessor>
void particle_strategy<eorder,accessor>::m2m( Mcoeff_t &target, const Mcoeff_t &source, point r ) const noexcept
{
    std::array<cmplx,eorder> transfer = I(-r);
    for ( uint k = 0; k < eorder; ++k )
    {
        for ( uint l = 0; l <= k; ++l )
        {
            target[k] += transfer[k-l]*source[l];
        }
    }
}

template <uint eorder, typename accessor> inline
void particle_strategy<eorder,accessor>::m2l( Lcoeff_t &target, const Mcoeff_t &source, point r ) const noexcept
{
    real fac {1};
    std::array<cmplx,eorder> transfer = O(r);
    for ( uint l = 0; l < eorder; ++l )
    {
        for ( uint k = 0; k < eorder - l; ++k )
        {
            target[l] += fac*transfer[l+k]*source[k];
        }
        fac = -fac;
    }
}

template <uint eorder, typename accessor>
void particle_strategy<eorder,accessor>::l2l( Lcoeff_t &target, const Mcoeff_t &source, point r ) const noexcept
{
    std::array<cmplx,eorder> transfer = I(r);
    for ( uint l = 0; l < eorder; ++l )
    {
        for ( uint k = 0; k < eorder - l; ++k )
        {
            target[l] += transfer[k]*source[k+l];
        }
    } 
}

template <uint eorder, typename accessor>
template <typename target_iterator>
void particle_strategy<eorder,accessor>::m2p( target_iterator *const begin, target_iterator *const end,
                                              const Mcoeff_t& M, point xc ) noexcept
{
    constexpr real  fac { 1./(2.*math::pi) };
    for ( target_iterator *i = begin; i != end; ++i )
    {
        target_iterator t = *i;
        const point r = get(*t) - xc;
        const cmplx z { r.x, r.y };
        Mcoeff_t transfer = dO(r);

        for ( uint k = 0; k < eorder; ++k )
        {
            (get.u(*t)).x -= fac*(transfer[k]*M[k]).imag();
            (get.u(*t)).y -= fac*(transfer[k]*M[k]).real();
        }
    }
}

template <uint eorder,typename accessor>
template <typename target_iterator>
void particle_strategy<eorder,accessor>::l2p( target_iterator *const begin, target_iterator *const end,
                                            const Lcoeff_t &L, point xc ) noexcept
{
    constexpr real  fac { 1./(2.*math::pi) };
    for ( target_iterator *i = begin; i != end; ++i )
    {
        target_iterator t = *i;
        const point r = get(*t) - xc;
        Mcoeff_t  transfer =  dI(r);

        for ( uint k = 0; k < eorder; ++k )
        {
            (get.u(*t)).x -= fac*(transfer[k]*L[k]).imag();
            (get.u(*t)).y -= fac*(transfer[k]*L[k]).real();
        }
    }
}

template <uint eorder,typename accessor>
template <typename target_iterator, typename source_iterator>
void particle_strategy<eorder,accessor>::p2p( target_iterator *const tbegin, target_iterator *const tend,
                                              source_iterator *const sbegin, source_iterator *const send ) noexcept
{
    constexpr real  fac { 1./(2.*math::pi) };
    for ( target_iterator *i = tbegin; i != tend; ++i )
    {
        target_iterator t = *i;
        for ( const source_iterator *j = sbegin; j != send; ++j )
        {
            source_iterator s = *j;

            const real rx = get(*t).x - get(*s).x;
            const real ry = get(*t).y - get(*s).y;
            const real  G = get.G(*s).z;
            const real r2 = rx*rx + ry*ry;

            if ( r2 == 0 ) continue;

            const real moll = -std::expm1( -r2/(epsilon*epsilon) );
            //const real reps = r2/(epsilon*epsilon);
            //const real moll = 1.-(1.-2.*reps+reps*reps/2.)*std::exp( -reps );
            (get.u(*t)).x -= fac*moll*G*ry / r2;
            (get.u(*t)).y += fac*moll*G*rx / r2;
        }
    }
}

template <uint eorder, typename accessor> inline
bool particle_strategy<eorder,accessor>::mac( point Acentre, real Aradius, point Bcentre, real Bradius ) const noexcept
{
    return (Aradius+Bradius)/length( Acentre - Bcentre ) < theta_max;
}

template <uint eorder,typename accessor>
template <typename iterator> inline
geometry::point particle_strategy<eorder,accessor>::pos( iterator it ) const noexcept
{
    return get(*it); 
}

template <uint eorder,typename accessor>
template <typename iterator> inline
void particle_strategy<eorder,accessor>::bounding_box( iterator it, geometry::point& min, geometry::point& max ) const noexcept
{
    min = max = get(*it);
}

template <uint eorder,typename accessor> inline
std::array<cmplx,eorder> particle_strategy<eorder,accessor>::I( point x ) const noexcept
{
    const cmplx z { x.x, x.y };
    std::array<cmplx,eorder> result;

    result[ 0 ] = { 1, 0 };
    for ( uint k = 1; k < eorder; ++k )
    {
        result[ k ] = result[ k - 1 ] * z / ((real) k);
    }
    return result;
}

template <uint eorder, typename accessor> inline
std::array<cmplx,eorder> particle_strategy<eorder,accessor>::dI( point x ) const noexcept
{
    const cmplx z { x.x, x.y };
    std::array<cmplx,eorder> result;

    result[ 0 ] = { 0, 0 };
    result[ 1 ] = { 1, 0 };
    for ( uint k = 2; k < eorder; ++k )
    {
        result[ k ] = result[ k - 1 ] * z / ((real) (k-1));
    }
    return result;
}

template <uint eorder, typename accessor> inline
std::array<cmplx,eorder> particle_strategy<eorder,accessor>::O( point x ) const noexcept
{
    const cmplx z   { x.x, x.y };
    std::array<cmplx,eorder> result;

    result[ 0 ] = -std::log( z );
    result[ 1 ] = 1./z;
    for ( uint k = 2; k < eorder; ++k )
    {
        result[k] = result[k-1] * ((real) (k-1)) / z;
    }
    return result;
}

template <uint eorder, typename accessor> inline
std::array<cmplx,eorder> particle_strategy<eorder,accessor>::dO( point x ) const noexcept
{
    const cmplx z   { x.x, x.y };
    std::array<cmplx,eorder> result;

    result[0] = -1./z;
    for ( uint k = 1; k < eorder; ++k )
    {
        result[k] = result[k-1] * ((real) k) / z;
    }
    return result;
}

}

