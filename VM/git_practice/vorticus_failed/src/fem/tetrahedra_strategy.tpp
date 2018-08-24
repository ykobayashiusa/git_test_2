/*
 * Copyright (C) 2016 Matthias Kirchhart
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

namespace fem
{

template <uint eorder>
template <typename source_iter>               
void tetrahedra_strategy<eorder>::p2m
(
    Mcoeff_t &M, point xc, const source_iter *const begin,
                           const source_iter *const end
) const 
{
    using math::dR;

    tetrahedral_quadrule qrule = get_tetrahedral_quadrule( eorder + 1 );

    for ( const source_iter *i = begin; i != end; ++i )
    {
        const tetrahedron<1>& t = **i;
        const real vol_elem = t.vol_elem( bary3d { .25, .25, .25, .25 } );

        for ( size_t q = 0; q < qrule.size(); ++q )
        {
            const bary3d pos = qrule[ q ].b;
            const real  w    = qrule[ q ].w * vol_elem / 6;

            const point yy  = t.Chi( pos );
            const point val = omega->eval( t, pos );

            const auto f   = dR<eorder>( yy - xc );

            for ( size_t k = 0; k < M.num_coeff(); ++k )
                M.data[k] += w*cross_prod( conj(f.data[k]), -val );
        }
    }
}

template <uint eorder>
template <typename source_iter>               
void tetrahedra_strategy<eorder>::p2l
(
    Mcoeff_t &L, point xc, const source_iter *const begin,
                           const source_iter *const end
) const
{
    using math::dS;

    tetrahedral_quadrule qrule = get_tetrahedral_quadrule( eorder + 1 );

    for ( const source_iter *i = begin; i != end; ++i )
    {
        const tetrahedron<1>& t = **i;
        const real vol_elem = t.vol_elem( bary3d { .25, .25, .25, .25 } );

        for ( size_t q = 0; q < qrule.size(); ++q )
        {
            const bary3d pos = qrule[ q ].b;
            const real  w    = qrule[ q ].w * vol_elem / 6;

            const point yy  = t.Chi( pos );
            const point val = omega->eval( t, pos );

            const auto f    = dS<eorder>( yy - xc );

            for ( size_t k = 0; k < L.num_coeff(); ++k )
                L.data[k] += w*cross_prod( f.data[k], -val );
        }
    }
}

template <uint eorder>
void tetrahedra_strategy<eorder>::m2m( Mcoeff_t &target, const Mcoeff_t &source, point r ) const noexcept
{
    using math::R;
    using math::idx;

    auto transfer = R<eorder>( -r );
    for ( int n = 0; n < (int) eorder; ++n )
    {
        for ( int m = 0; m <= n; ++m )
        {
            for ( int nd = 0; nd <= n; ++nd )
            {
                int min_md = std::max(-nd, m - (n-nd));
                int max_md = std::min( nd, m + (n-nd));
                for ( int md = min_md; md <= max_md; ++md )
                {
                    target.data[ idx(n,m) ] += source( nd, md )*conj(transfer(n-nd,m-md));
                }
            }
        }
    }
}

template <uint eorder>
void tetrahedra_strategy<eorder>::m2l( Lcoeff_t &target, const Mcoeff_t &source, point r ) const noexcept
{
    using math::idx;
    using math::S;
 
    const auto transfer = S<eorder>( r );

    for ( int n = 0; n < (int) eorder; ++n )
    {
        for ( int m = 0; m <= n; ++m )
        {
            cmplx_point tmp;
            for ( int nd = 0; nd < (int) eorder - n; ++nd )
            {
                for ( int md = -nd; md <= nd; ++md )
                {
                    tmp += source(nd,md)*transfer(n+nd,m+md);
                }
            }
            if ( n % 2 ) target.data[ idx(n,m) ] -= tmp;
            else         target.data[ idx(n,m) ] += tmp;
        }
    }
}

template <uint eorder>
void tetrahedra_strategy<eorder>::l2l( Lcoeff_t &target, const Mcoeff_t &source, point r ) const noexcept
{
    using math::idx;
    using math::R;
    const auto transfer = R<eorder>( r );

    for ( int n = 0; n < (int) eorder; ++n )
    {
        for ( int m = 0; m <= n; ++m )
        {
            for ( int nd = n; nd < (int) eorder; ++nd )
            {
                int min_md = std::max(-nd, m - (nd-n));
                int max_md = std::min( nd, m + (nd-n));
                for ( int md = min_md; md <= max_md; ++md )
                {
                    target.data[ idx(n,m) ] += source( nd,md )*conj(transfer(nd-n,md-m));
                }
            }
        }
    }
}

template <uint eorder>
template <typename target_iter>
void tetrahedra_strategy<eorder>::m2p( target_iter *const begin, target_iter *const end,
                                       const Mcoeff_t& M, point xc ) const
{
    using math::S;
    using math::pifac;

    for ( target_iter *it = begin; it != end; ++it )
    {
        const point r = (*it)->x - xc;
        const auto transfer = S<eorder>(r);

        point result {};
        for ( int n = 0; n < (int) eorder; ++n )
        {
            for ( int m = -n; m <= n; ++m )
            {
                result += (M(n,m)*transfer(n,m)).real()*pifac;
            }
        }
        (*it)->u += result;
    } 
}

template <uint eorder>
template <typename target_iter>
void tetrahedra_strategy<eorder>::l2p( target_iter *const begin, target_iter *const end,
                                       const Lcoeff_t &L, point xc ) const
{
    using math::R;
    using math::pifac;

    for ( target_iter *it = begin; it != end; ++it )
    {
        const point r = (*it)->x - xc;
        const auto transfer = R<eorder>(r);

        point result {};
        for ( int n = 0; n < (int) eorder; ++n )
        {
            for ( int m = -n; m <= n; ++m )
            {
                result += (L(n,m)*conj(transfer(n,m))).real()*pifac;
            }
        }
        (*it)-> u += result;
    } }

template <uint eorder>
template <typename target_iter, typename source_iter>
void tetrahedra_strategy<eorder>::p2p(       target_iter *const tbegin,       target_iter *const tend,
                                       const source_iter *const sbegin, const source_iter *const send ) const
{
    for ( const source_iter *j = sbegin; j != send; ++j )
    {
        const tetrahedron<1> &t = **j;
        const std::array<point,4> X { t.get_node(0), t.get_node(1),
                                      t.get_node(2), t.get_node(3) };
        const std::array<point,4> O { omega->eval(X[0]), omega->eval(X[1]),
                                      omega->eval(X[2]), omega->eval(X[3]) };

        for ( target_iter *i = tbegin; i != tend; ++i )
        {
            (*i)->u += biot_savart<1>( X, O, (*i)->x );
        }
    }
}

template <uint eorder>
bool tetrahedra_strategy<eorder>::mac( point Acentre, real Aradius, point Bcentre, real Bradius ) const noexcept
{
    return (Aradius+Bradius)/length( Acentre - Bcentre ) < 0.8;
}

template <uint eorder>
template <typename source_iter>
point tetrahedra_strategy<eorder>::pos( source_iter it ) const noexcept
{
    return it->pos();
}


template <uint eorder>
template <typename source_iter>
void tetrahedra_strategy<eorder>::bounding_box( source_iter it, geometry::point& min, geometry::point& max ) const noexcept
{
    it->bounding_box( min, max );
}

}

