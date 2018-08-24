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

namespace bem
{

template <uint eorder, uint gorder, uint aorder>
panel_strategy<eorder,gorder,aorder>::panel_strategy( const grid_t& g ):
elem_matrices ( g.size() ), partial_results ( g.size() )
{
    std::vector<const triangle<gorder>*> trias;

    // Create all entries.
    for ( const triangle<gorder>& t: g )
    {
        elem_matrices  [ t.pos() ];
        partial_results[ t.id()  ];
        trias.push_back( &t );
    }

    // Compute the self-interactions.
    #pragma omp parallel for
    for ( size_t i = 0; i < trias.size(); ++i )
    {
        const triangle<gorder>& t = *trias[ i ];
        get_elem_mat( t, t );
    }
}

template <uint eorder, uint gorder, uint aorder>
void panel_strategy<eorder,gorder,aorder>::apply( const grid_t& grid,
                                                  const grid_function<real,gorder,aorder>& source,
                                                        grid_function<real,gorder,aorder>& target )
{
    source_strengths = &source;

    for ( auto& value_pair: target )
    {
        value_pair.second = 0;
    }

    for ( auto& tmp: partial_results )
    {
        auto& vals = tmp.second;
        for ( size_t k = 0; k < vals.size(); ++k )
        {
            vals[ k ] = 0;
        }
    }

    fmm::fmm( *this, grid.begin(), grid.end(), grid.begin(), grid.end() );
    
    for ( const auto& tmp: partial_results )
    {
        const auto& vals          = tmp.second;
        const triangle<gorder>& t = grid.get_multi_grid().get_triangle( tmp.first );
        const auto& dof_map       = target.get_dofs( t );

        for ( size_t k = 0; k < vals.size(); ++k )
        {
            target( dof_map[ k ] ) += vals[ k ];
        }
    }
}

template <uint eorder, uint gorder, uint aorder>
void panel_strategy<eorder,gorder,aorder>::p2m
(
    Mcoeff_t &M, point xc, const source_iter *const begin,
                           const source_iter *const end
) const noexcept
{
    using math::R;
    using shapefcts2d::N;
    using shapefcts2d::apply;
    using geometry::cmplx_point;


    const uint  quad_order   = eorder + aorder + 1;
    const uint  quad_num     = dunavant_num_points[ quad_order ];
    const bary *quad_points  = dunavant_points [ quad_order ];
    const real *quad_weights = dunavant_weights[ quad_order ];

    for ( const source_iter *i = begin; i != end; ++i )
    {
        const triangle<gorder>& t = **i;
        const auto dof_map = source_strengths->get_dofs( t );
 
        std::array<real,dof_map.size()> dofs;
        for ( size_t j = 0; j < dofs.size(); ++j )
        {
            dofs[ j ] = (*source_strengths)( dof_map[j] );
        }

        for ( uint q = 0; q < quad_num; ++q )
        {
            const bary pos = quad_points [ q ];
            const real w   = quad_weights[ q ];

            const point yy = t.Chi( pos );
            const real  gy = t.surf_elem( pos );

            const auto phi = rotN<gorder,aorder>( t, pos );
            const auto f   = R<eorder>( yy - xc );

            const auto val = apply<aorder,point>( phi, dofs );

            for ( size_t k = 0; k < M.num_coeff(); ++k )
                M.data[k] += 0.5*w*gy*conj(f.data[k])*val;
        }
    }
}

template <uint eorder, uint gorder, uint aorder>
void panel_strategy<eorder,gorder,aorder>::p2l
(
    Mcoeff_t &L, point xc, const source_iter *const begin,
                           const source_iter *const end
) const noexcept
{
    using math::S;
    using shapefcts2d::N;
    using shapefcts2d::apply;
    using geometry::cmplx_point;


    const uint  quad_order   = eorder + aorder + 1;
    const uint  quad_num     = dunavant_num_points[ quad_order ];
    const bary *quad_points  = dunavant_points [ quad_order ];
    const real *quad_weights = dunavant_weights[ quad_order ];

    for ( const source_iter *i = begin; i != end; ++i )
    {
        const triangle<gorder>& t = **i;
        const auto dof_map = source_strengths->get_dofs( t );
 
        std::array<real,dof_map.size()> dofs;
        for ( size_t j = 0; j < dofs.size(); ++j )
        {
            dofs[ j ] = (*source_strengths)( dof_map[j] );
        }

        for ( uint q = 0; q < quad_num; ++q )
        {
            const bary pos = quad_points [ q ];
            const real w   = quad_weights[ q ];

            const point yy = t.Chi( pos );
            const real  gy = t.surf_elem( pos );

            const auto phi = rotN<gorder,aorder>( t, pos );
            const auto f   = S<eorder>( yy - xc );

            const auto val = apply<aorder,point>( phi, dofs );

            for ( size_t k = 0; k < L.num_coeff(); ++k )
                L.data[k] += 0.5*w*gy*f.data[k]*val;
        }
    }
}

template <uint eorder, uint gorder, uint aorder>
void panel_strategy<eorder,gorder,aorder>::m2m( Mcoeff_t &target, const Mcoeff_t &source, point r ) const noexcept
{
    using math::R;
    using math::idx;

    auto transfer = R<eorder>( -r );

    constexpr int order = eorder;
    for ( int n = 0; n < order; ++n )
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

template <uint eorder, uint gorder, uint aorder>
void panel_strategy<eorder,gorder,aorder>::m2l( Lcoeff_t &target, const Mcoeff_t &source, point r ) const noexcept
{
     using math::idx;
     using math::S;
 
    const auto transfer = S<eorder>( r );
    constexpr int order = eorder;
    for ( int n = 0; n < order; ++n )
    {
        for ( int m = 0; m <= n; ++m )
        {
            cmplx_point tmp;
            for ( int nd = 0; nd < order - n; ++nd )
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

/*
template <uint eorder, uint gorder, uint aorder>
void panel_strategy<eorder,gorder,aorder>::m2l( Lcoeff_t &target, const Mcoeff_t &source, point r ) const noexcept
{
    using math::idx;
    using math::S;

    constexpr int order = eorder;
    constexpr int o     = eorder;
    const real d = r.r()/order;
    const auto transfer = S<eorder>( r/d );

    arma::cx_mat::fixed<2*eorder,2*eorder> fft_x, 
                                           fft_y,
                                           fft_z,
                                           transfer_fft;
    fft_x.fill(0);
    fft_y.fill(0);
    fft_z.fill(0);
    transfer_fft.fill(0);

    real dn = 1;
    for ( int n = 0; n < order; ++n )
    {
        for ( int m = -n; m <= n; ++m )
        {
            fft_x( o+n, o+m ) = source( n, m ).x / dn;
            fft_y( o+n, o+m ) = source( n, m ).y / dn;
            fft_z( o+n, o+m ) = source( n, m ).z / dn;

            transfer_fft( o-n, o-m ) = transfer( n, m );
        }
        dn *= d;
    }

    fftw_execute_dft( fwd_plan, reinterpret_cast<fftw_complex*>(fft_x.memptr()), reinterpret_cast<fftw_complex*>(fft_x.memptr()) );
    fftw_execute_dft( fwd_plan, reinterpret_cast<fftw_complex*>(fft_y.memptr()), reinterpret_cast<fftw_complex*>(fft_y.memptr()) );
    fftw_execute_dft( fwd_plan, reinterpret_cast<fftw_complex*>(fft_z.memptr()), reinterpret_cast<fftw_complex*>(fft_z.memptr()) );
    fftw_execute_dft( fwd_plan, reinterpret_cast<fftw_complex*>(transfer_fft.memptr()), reinterpret_cast<fftw_complex*>(transfer_fft.memptr()) );

    fft_x %= transfer_fft;
    fft_y %= transfer_fft;
    fft_z %= transfer_fft;

    fftw_execute_dft( bwd_plan, reinterpret_cast<fftw_complex*>(fft_x.memptr()), reinterpret_cast<fftw_complex*>(fft_x.memptr()) );
    fftw_execute_dft( bwd_plan, reinterpret_cast<fftw_complex*>(fft_y.memptr()), reinterpret_cast<fftw_complex*>(fft_y.memptr()) );
    fftw_execute_dft( bwd_plan, reinterpret_cast<fftw_complex*>(fft_z.memptr()), reinterpret_cast<fftw_complex*>(fft_z.memptr()) );

    constexpr real fac = 1./(4*order*order);
    target.data[ 0 ].x += (fac/d)*fft_x(0,0);
    target.data[ 0 ].y += (fac/d)*fft_y(0,0);
    target.data[ 0 ].z += (fac/d)*fft_z(0,0);

    dn = d*d;
    for ( int n = 1; n < order; ++n )
    {
        for ( int m = 0; m <= n; ++m )
        {
            real sign = (n+m) % 2 ? -1:1;
            target.data[ idx(n,m) ].x += (fac/dn)*sign*conj(fft_x( 2*order - n, m ));
            target.data[ idx(n,m) ].y += (fac/dn)*sign*conj(fft_y( 2*order - n, m ));
            target.data[ idx(n,m) ].z += (fac/dn)*sign*conj(fft_z( 2*order - n, m ));
        }
        dn *= d;
    }
}
*/

template <uint eorder, uint gorder, uint aorder>
void panel_strategy<eorder,gorder,aorder>::l2l( Lcoeff_t &target, const Mcoeff_t &source, point r ) const noexcept
{
    using math::idx;
    using math::R;
    const auto transfer = R<eorder>( r );

    constexpr int order = eorder;
    for ( int n = 0; n < order; ++n )
    {
        for ( int m = 0; m <= n; ++m )
        {
            for ( int nd = n; nd < order; ++nd )
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

template <uint eorder, uint gorder, uint aorder>
void panel_strategy<eorder,gorder,aorder>::m2p( target_iter *const begin, target_iter *const end,
                                                const Mcoeff_t& M, point xc ) noexcept
{
    const uint  quad_order   = eorder + aorder + 1;
    const uint  quad_num     = dunavant_num_points[ quad_order ];
    const bary *quad_points  = dunavant_points [ quad_order ];
    const real *quad_weights = dunavant_weights[ quad_order ];
   
    for ( target_iter *it = begin; it != end; ++it )
    {
        const triangle<gorder>& t = **it;
        shapefcts2d::vals<aorder>& target = partial_results[ t.id() ];

        for ( uint q = 0; q < quad_num; ++q )
        {
            const bary pos = quad_points [ q ];
            const real w   = quad_weights[ q ];

            const point x = t.Chi      ( pos );
            const real gx = t.surf_elem( pos );

            const auto  rot  = rotN<gorder,aorder>( t, pos );
            const point val  = eval_m_exp( M, x - xc );

            for ( size_t i = 0; i < target.size(); ++i )
            {
                target[i] += 0.5*w*gx*scal_prod(rot[i],val);
            }
        }
    } 
}

template <uint eorder, uint gorder, uint aorder>
void panel_strategy<eorder,gorder,aorder>::l2p( target_iter *const begin, target_iter *const end,
                                                const Lcoeff_t &L, point xc ) noexcept
{
    const uint  quad_order   = eorder + aorder + 1;
    //const uint  quad_order   = 3; // Outer integration only needs to be of exactness degree 3!
    const uint  quad_num     = dunavant_num_points[ quad_order ];
    const bary *quad_points  = dunavant_points [ quad_order ];
    const real *quad_weights = dunavant_weights[ quad_order ];
   
    for ( target_iter *it = begin; it != end; ++it )
    {
        const triangle<gorder>& t = **it;
        shapefcts2d::vals<aorder>& target = partial_results[ t.id() ];

        for ( uint q = 0; q < quad_num; ++q )
        {
            const bary pos = quad_points [ q ];
            const real w   = quad_weights[ q ];

            const point x = t.Chi      ( pos );
            const real gx = t.surf_elem( pos );

            const auto  rot  = rotN<gorder,aorder>( t, pos );
            const point val  = eval_l_exp( L, x - xc );

            for ( size_t i = 0; i < target.size(); ++i )
                target[i] += 0.5*w*gx*scal_prod(rot[i],val);
        }
    } 
}

template <uint eorder, uint gorder, uint aorder> inline
auto panel_strategy<eorder,gorder,aorder>::eval_l_exp( const Lcoeff_t &L, point r )
const noexcept -> point
{
    using math::R;
    using math::pifac;

    point result {};
    const auto transfer = R<eorder>(r);

    constexpr int order = eorder;
    for ( int n = 0; n < order; ++n )
    {
        for ( int m = -n; m <= n; ++m )
        {
            result += (L(n,m)*conj(transfer(n,m))).real()*pifac;
        }
    }

    return result;
}

template <uint eorder, uint gorder, uint aorder> inline
auto panel_strategy<eorder,gorder,aorder>::eval_m_exp( const Mcoeff_t &M, point r )
const noexcept -> point
{
    using math::S;
    using math::pifac;

    point result {};
    const auto transfer = S<eorder>(r);

    constexpr int order = eorder;
    for ( int n = 0; n < order; ++n )
    {
        for ( int m = -n; m <= n; ++m )
        {
            result += (M(n,m)*transfer(n,m)).real()*pifac;
        }
    }

    return result;
}

template <uint eorder, uint gorder, uint aorder>
void panel_strategy<eorder,gorder,aorder>::p2p(       target_iter *const tbegin,       target_iter *const tend,
                                                const source_iter *const sbegin, const source_iter *const send ) noexcept
{
    arma::vec::fixed< shapefcts2d::num<aorder>() > rhs;

    for ( target_iter *i = tbegin; i != tend; ++i )
    {
        const triangle<gorder> &t1 = **i;
        const auto  dofs1 = source_strengths->get_dofs( t1 );

        for ( const source_iter *j = sbegin; j != send; ++j )
        {
            const triangle<gorder> &t2 = **j;
            const auto  dofs2 = source_strengths->get_dofs( t2 );

            const elem_mat_t& elem_mat = get_elem_mat( t1, t2 );

            for ( size_t k = 0; k < dofs2.size(); ++k )
            {
                rhs( k ) = (*source_strengths)( dofs2[ k ] );
            }

            rhs = elem_mat*rhs;
        
            shapefcts2d::vals<aorder>& target = partial_results[ t1.id() ];
            for ( size_t k = 0; k < dofs1.size(); ++k )
            {
                target[ k ] += rhs(k);
            }
        }
    }
}

template <uint eorder, uint gorder, uint aorder>
bool panel_strategy<eorder,gorder,aorder>::mac( point Acentre, real Aradius, point Bcentre, real Bradius ) const noexcept
{
    return (Aradius+Bradius)/length( Acentre - Bcentre ) < 0.8;
}

template <uint eorder, uint gorder, uint aorder>
auto panel_strategy<eorder,gorder,aorder>::get_elem_mat( const triangle<gorder> &t1,
                                                         const triangle<gorder> &t2 )
-> const elem_mat_t&
{
    if ( elem_matrices[ t1.pos() ].find( t2.pos() ) == elem_matrices[ t1.pos() ].end() )
    {
        elem_matrices[ t1.pos() ][ t2.pos() ] = double_integral( W_kernel<gorder,aorder>, t1, t2 );
    }

    return elem_matrices[ t1.pos() ][ t2.pos() ];
}

}

