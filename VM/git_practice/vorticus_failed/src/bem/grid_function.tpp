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


template <typename T, uint gorder, uint aorder>
grid_function<T,gorder,aorder>::grid_function( const grid<gorder> &g ):
level_ { g.get_level() }
{
    using std::make_pair;

    coords_.reserve( g.size() );
    const T default_value {};
    for ( const auto &t: g )
    {
        coords_t& map = coords_[ t.id() ];

        // Corner dofs.
        values_.insert( make_pair(t.get_node(0), default_value) ); map[ 0 ] = t.get_node(0);
        values_.insert( make_pair(t.get_node(1), default_value) ); map[ 1 ] = t.get_node(1);
        values_.insert( make_pair(t.get_node(2), default_value) ); map[ 2 ] = t.get_node(2);

        // Edge dofs.
        for ( unsigned char e_num { 0 }; e_num < 3; ++e_num )
        {
            const edge<gorder>& e { t.get_edge(e_num) };
            
            size_t offset;
            switch ( e_num )
            {
            case 0: offset =   aorder + 1; break;
            case 1: offset = 2*aorder + 0; break;
            case 2: offset =            2; break;
            } 

            bool reverse;
            if ( (e.has_triangle(0) && &e.get_triangle(0) == &t) ||
                 (e.has_triangle(2) && &e.get_triangle(2) == &t) )
            {
                reverse = false;
            }
            else
            {
                reverse = true;
            }
            
            for ( size_t k { 1 }; k < aorder; ++k )
            {
                const real pos { ((real) k)/((real) aorder) };
                values_.insert( make_pair(e.Chi(pos), default_value) );

                const size_t tria_dof_no = reverse ? offset+(aorder-k) : offset+k;
                map[ tria_dof_no ] = e.Chi(pos);
            }
        }

        // Inner dofs.
        constexpr auto pos = lattice::lattice<aorder>();
        for ( size_t k { 3*aorder }; k < lattice::size<aorder>(); ++k )
        {
            values_.insert( make_pair(t.Chi(pos[k]), default_value) );
            map[ k ] = t.Chi(pos[k]);
        }
    }
}

template <typename T, uint gorder, uint aorder>
T grid_function<T,gorder,aorder>::eval( const triangle<gorder>& t, bary pos ) const
{
    using namespace shapefcts2d;
    coeffs<aorder,T> coeffs;
    const coords_t& map = coords_.find( t.id() )->second;

    for ( size_t i { 0 }; i < map.size(); ++i )
    {
        coeffs[ i ] = values_.find( map[ i ] )->second;
    }
   
    return apply<aorder,T>( coeffs, N<aorder>(pos) );
}

template <typename T, uint gorder, uint aorder> inline
T grid_function<T,gorder,aorder>::eval( point pos ) const
{
    return values_.find( pos )->second;
}

template <typename T, uint gorder, uint aorder> inline
T& grid_function<T,gorder,aorder>::operator()( point pos ) noexcept
{
    return values_.find( pos )->second;
}

template <typename T, uint gorder, uint aorder> inline
const T& grid_function<T,gorder,aorder>::operator()( point pos ) const noexcept
{
    return values_.find( pos )->second;
}

template <typename T, uint gorder, uint aorder> inline
auto grid_function<T,gorder,aorder>::begin() noexcept -> iterator
{
    return values_.begin();
}

template <typename T, uint gorder, uint aorder> inline
auto grid_function<T,gorder,aorder>::begin() const noexcept -> const_iterator
{
    return values_.begin();
}

template <typename T, uint gorder, uint aorder> inline
auto grid_function<T,gorder,aorder>::end() noexcept -> iterator
{
    return values_.end();
}

template <typename T, uint gorder, uint aorder> inline
auto grid_function<T,gorder,aorder>::end() const noexcept -> const_iterator
{
    return values_.end();
}

template <typename T, uint gorder, uint aorder> inline
auto grid_function<T,gorder,aorder>::size() const noexcept -> size_type
{
    return values_.size();
}

template <typename T, uint gorder, uint aorder> inline
auto grid_function<T,gorder,aorder>::get_dofs( const triangle<gorder>& t ) const noexcept -> coords_t
{
    return coords_.find( t.id() )->second;
}

template <typename T, uint gorder, uint aorder> inline
void grid_function<T,gorder,aorder>::operator+=( const grid_function& rhs )
{
    for ( auto& val: values_ )
        val.second += rhs( val.first );
}

template <typename T, uint gorder, uint aorder> inline
void grid_function<T,gorder,aorder>::operator-=( const grid_function& rhs )
{
    for ( auto& val: values_ )
        val.second -= rhs( val.first );
}

template <typename T, uint gorder, uint aorder> inline
void grid_function<T,gorder,aorder>::operator*=( const real rhs )
{
    for ( auto& val: values_ )
        val.second *= rhs;
}

template <typename T, uint gorder, uint aorder> inline
void grid_function<T,gorder,aorder>::operator/=( const real rhs )
{
    for ( auto& val: values_ )
        val.second /= rhs;
}

template <typename T, uint gorder, uint aorder> inline
auto grid_function<T,gorder,aorder>::operator=( const T& rhs ) -> grid_function&
{
    for ( auto& val: values_ )
        val.second = rhs;
    return *this;
}

template <uint gorder, uint aorder>
prolongation<gorder,aorder>::prolongation( const typename multigrid<gorder>::grid& gcoarse,
                                           const typename multigrid<gorder>::grid& gfine )
{
    grid_function<real,gorder,aorder> coarse( gcoarse );
    grid_function<real,gorder,aorder> fine  ( gfine   );

    for ( auto& t: gcoarse )
    {
        treat_triangle(t,coarse,fine);
    }
}

template <uint gorder, uint aorder>
template <typename T>
void prolongation<gorder,aorder>::prolongate( const grid_function<T,gorder,aorder>& coarse,
                                                    grid_function<T,gorder,aorder>&   fine ) const
{
    for ( auto& val: fine.values_ )
        val.second = 0;

    for ( const auto& entry: matrix_ )
    {
        const point coarse_node { entry.first.first  };
        const point   fine_node { entry.first.second };
        const real        value { entry.second };

        fine( fine_node ) += coarse( coarse_node )*value;
    }
}

template <uint gorder, uint aorder>
template <typename T>
void prolongation<gorder,aorder>::restrict(       grid_function<T,gorder,aorder>& coarse,
                                            const grid_function<T,gorder,aorder>&   fine ) const
{
    for ( auto& val: coarse.values_ )
        val.second = 0;

    for ( const auto& entry: matrix_ )
    {
        const point coarse_node { entry.first.first  };
        const point   fine_node { entry.first.second };
        const real        value { entry.second };

        coarse( coarse_node ) += fine( fine_node )*value;
    }
}

template <uint gorder, uint aorder>
void prolongation<gorder,aorder>::treat_triangle( const triangle<gorder>& t,
                                                  const scalar& coarse, const scalar& fine )
{
    const auto  ref_status = t.get_ref_status();
    
    switch ( ref_status.count() )
    {
    case 0: treat_no_ref          ( t, coarse, fine ); break;
    case 1: treat_single_irreg_ref( t, coarse, fine ); break;
    case 2: treat_double_irreg_ref( t, coarse, fine ); break;
    case 3: treat_regular_ref     ( t, coarse, fine ); break;
    }
}

template <uint gorder, uint aorder>
void prolongation<gorder,aorder>::treat_no_ref( const triangle<gorder>& t,
                                                const scalar &coarse, const scalar& )
{
    const auto& coarse_map = coarse.coords_.find( t.id() )->second;

    // Simply copy the values.
    for ( const point pos: coarse_map )
    {
        const idx_t idx { pos, pos };
        matrix_[ idx ] = 1;
    } 
}

template <uint gorder, uint aorder>
void prolongation<gorder,aorder>::treat_single_irreg_ref( const triangle<gorder>& t,
                                                          const scalar &coarse, const scalar &fine )
{
    const auto& coarse_map   = coarse.coords_.find( t.id() )->second;
    const auto  ref_status   = t.get_ref_status();
    const auto&   fine_map_0 = fine.coords_.find( t.get_child_id( 0 ) )->second;
    const auto&   fine_map_1 = fine.coords_.find( t.get_child_id( 1 ) )->second;

    bary t0_0, t0_1, t0_2;
    bary t1_0, t1_1, t1_2;
    if ( ref_status.test(0) )
    {
        t0_0 = bary {  1, .0, .0 }; t0_1 = bary { .0,  1, .0 }; t0_2 = bary { .0, .5, .5 };
        t1_0 = bary { .0, .0,  1 }; t1_1 = bary {  1, .0, .0 }; t1_2 = bary { .0, .5, .5 };
    }
    else if ( ref_status.test(1) )
    {
        t0_0 = bary {  1, .0, .0 }; t0_1 = bary { .0,  1, .0 }; t0_2 = bary { .5, .0, .5 };
        t1_0 = bary { .0,  1, .0 }; t1_1 = bary { .0, .0,  1 }; t1_2 = bary { .5, .0, .5 };
    }
    else // ( ref_status.test(2) )
    {
        t0_0 = bary { .0, .0,  1 }; t0_1 = bary {  1, .0, .0 }; t0_2 = bary { .5, .5, .0 };
        t1_0 = bary { .0,  1, .0 }; t1_1 = bary { .0, .0,  1 }; t1_2 = bary { .5, .5, .0 };
    }
    
    constexpr auto positions = lattice::lattice<aorder>();
    for ( size_t k { 0 }; k < positions.size(); ++k )
    {
        const bary pos = positions[ k ];
        const bary pos0 { pos.z0*t0_0 + pos.z1*t0_1 + pos.z2*t0_2 };
        const bary pos1 { pos.z0*t1_0 + pos.z1*t1_1 + pos.z2*t1_2 };

        using shapefcts2d::N;
        const auto vals0 = N<aorder>( pos0 );
        const auto vals1 = N<aorder>( pos1 );

        for ( size_t j { 0 }; j < positions.size(); ++j )
        {
            const idx_t idx0 { coarse_map[ j ], fine_map_0[ k ] };
            const idx_t idx1 { coarse_map[ j ], fine_map_1[ k ] };

            matrix_[ idx0 ] = vals0[ j ];
            matrix_[ idx1 ] = vals1[ j ];
        }
    }
}

template <uint gorder, uint aorder>
void prolongation<gorder,aorder>::treat_double_irreg_ref( const triangle<gorder>& t,
                                                          const scalar &coarse, const scalar &fine )
{
    const auto& coarse_map   = coarse.coords_.find( t.id() )->second;
    const auto  ref_status   = t.get_ref_status();
    const auto&   fine_map_0 =   fine.coords_.find( t.get_child_id( 0 ) )->second;
    const auto&   fine_map_c =   fine.coords_.find( t.get_child_id( 1 ) )->second;
    const auto&   fine_map_1 =   fine.coords_.find( t.get_child_id( 2 ) )->second;

    bary t0_0, t0_1, t0_2;
    bary tc_0, tc_1, tc_2;
    bary t1_0, t1_1, t1_2;

    if ( ! ref_status.test(2) ) // 011
    {
        t0_0 = bary {  1, .0, .0 }; t0_1 = bary { .0,  1, .0 }; t0_2 = bary { .0, .5, .5 };
        tc_0 = bary {  1, .0, .0 }; tc_1 = bary { .0, .5, .5 }; tc_2 = bary { .5, .0, .5 };
        t1_0 = bary { .0, .0,  1 }; t1_1 = bary { .5, .0, .5 }; t1_2 = bary { .0, .5, .5 };
    }
    else if ( ! ref_status.test(1) ) // 101
    {
        t0_0 = bary { .0,  1, .0 }; t0_1 = bary { .0, .5, .5 }; t0_2 = bary { .5, .5, .0 };
        tc_0 = bary {  1, .0, .0 }; tc_1 = bary { .5, .5, .0 }; tc_2 = bary { .0, .5, .5 };
        t1_0 = bary { .0, .0,  1 }; t1_1 = bary {  1, .0, .0 }; t1_2 = bary { .0, .5, .5 };
    }
    else // 110
    {
        t0_0 = bary {  1, .0, .0 }; t0_1 = bary { .5, .5, .0 }; t0_2 = bary { .5, .0, .5 };
        tc_0 = bary { .0,  1, .0 }; tc_1 = bary { .5, .0, .5 }; tc_2 = bary { .5, .5, .0 };
        t1_0 = bary { .0, .0,  1 }; t1_1 = bary { .5, .0, .5 }; t1_2 = bary { .0,  1, .0 };
    }
    
    constexpr auto positions = lattice::lattice<aorder>();
    for ( size_t k { 0 }; k < positions.size(); ++k )
    {
        const bary pos = positions[ k ];
        const bary pos0 { pos.z0*t0_0 + pos.z1*t0_1 + pos.z2*t0_2 };
        const bary posc { pos.z0*tc_0 + pos.z1*tc_1 + pos.z2*tc_2 };
        const bary pos1 { pos.z0*t1_0 + pos.z1*t1_1 + pos.z2*t1_2 };

        using shapefcts2d::N;
        const auto vals0 = N<aorder>( pos0 );
        const auto valsc = N<aorder>( posc );
        const auto vals1 = N<aorder>( pos1 );

        for ( size_t j { 0 }; j < positions.size(); ++j )
        {
            const idx_t idx0 { coarse_map[ j ], fine_map_0[ k ] };
            const idx_t idxc { coarse_map[ j ], fine_map_c[ k ] };
            const idx_t idx1 { coarse_map[ j ], fine_map_1[ k ] };

            matrix_[ idx0 ] = vals0[ j ];
            matrix_[ idxc ] = valsc[ j ];
            matrix_[ idx1 ] = vals1[ j ];
        }
    }
}

template <uint gorder, uint aorder>
void prolongation<gorder,aorder>::treat_regular_ref( const triangle<gorder>& t,
                                                     const scalar& coarse, const scalar& fine )
{
    const auto& coarse_map   = coarse.coords_.find( t.id() )->second;
    const auto&   fine_map_0 =   fine.coords_.find( t.get_child_id( 0 ) )->second;
    const auto&   fine_map_1 =   fine.coords_.find( t.get_child_id( 1 ) )->second;
    const auto&   fine_map_2 =   fine.coords_.find( t.get_child_id( 2 ) )->second;
    const auto&   fine_map_3 =   fine.coords_.find( t.get_child_id( 3 ) )->second;

    constexpr bary t0_0 { 1., 0., 0. }, t0_1 { .5, .5, .0 }, t0_2 { .5, .0, .5 };
    constexpr bary t1_0 { 0., 1., 0. }, t1_1 { .0, .5, .5 }, t1_2 { .5, .5, .0 };
    constexpr bary t2_0 { 0., 0., 1. }, t2_1 { .5, .0, .5 }, t2_2 { .0, .5, .5 };
    constexpr bary t3_0 { 0., .5, .5 }, t3_1 { .5, .0, .5 }, t3_2 { .5, .5, .0 };
    constexpr auto positions = lattice::lattice<aorder>();

    for ( size_t k { 0 }; k < positions.size(); ++k )
    {
        const bary pos = positions[ k ];
        const bary pos0 { pos.z0*t0_0 + pos.z1*t0_1 + pos.z2*t0_2 };
        const bary pos1 { pos.z0*t1_0 + pos.z1*t1_1 + pos.z2*t1_2 };
        const bary pos2 { pos.z0*t2_0 + pos.z1*t2_1 + pos.z2*t2_2 };
        const bary pos3 { pos.z0*t3_0 + pos.z1*t3_1 + pos.z2*t3_2 };

        using shapefcts2d::N;
        const auto vals0 = N<aorder>( pos0 );
        const auto vals1 = N<aorder>( pos1 );
        const auto vals2 = N<aorder>( pos2 );
        const auto vals3 = N<aorder>( pos3 );

        for ( size_t j { 0 }; j < positions.size(); ++j )
        {
            const idx_t idx0 { coarse_map[ j ], fine_map_0[ k ] };
            const idx_t idx1 { coarse_map[ j ], fine_map_1[ k ] };
            const idx_t idx2 { coarse_map[ j ], fine_map_2[ k ] };
            const idx_t idx3 { coarse_map[ j ], fine_map_3[ k ] };

            matrix_[ idx0 ] = vals0[ j ];
            matrix_[ idx1 ] = vals1[ j ];
            matrix_[ idx2 ] = vals2[ j ];
            matrix_[ idx3 ] = vals3[ j ];
        }
    }
}

}

