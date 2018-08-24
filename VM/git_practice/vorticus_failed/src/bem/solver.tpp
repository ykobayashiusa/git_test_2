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

template <uint eorder, uint gorder, uint aorder>
solver<eorder,gorder,aorder>::solver( multigrid<gorder>& mg ):
mg_ (mg) 
{
    for ( size_t l = 0; l <= mg.last_level(); ++l )
    {
             grids_.push_back( grid(mg,l) );
        strategies_.push_back( strategy(grids_[l]) );
           defects_.push_back( fct( grids_[l] ) );

        if ( l >= 1 )
             prols_.push_back( prol( grids_[ l - 1 ], grids_[ l ] ) );

        numberings_.push_back( std::unordered_map<point,size_t>() );
        auto& numb = numberings_.back();
        size_t counter { 0 };
        for ( const auto& p: defects_[ l ] ) numb[ p.first ] = counter++;

        diagonals_.push_back( math::vector( numb.size(), arma::fill::zeros ) );
        auto& diag = diagonals_[ l ];

        for ( const auto& t: grids_[l] )
        {
            const auto& elem_mat = strategies_[l].get_elem_mat( t, t );
            const auto& dofs = defects_[l].get_dofs( t );

            for ( size_t i = 0; i < elem_mat.n_rows; ++i )
                diag( numb[ dofs[ i ] ] ) += elem_mat( i, i );
        }
    }
}

template <uint eorder, uint gorder, uint aorder>
auto solver<eorder,gorder,aorder>::compute_rhs( const fct& g_N ) -> fct
{
    unsigned char level = g_N.get_level();
    math::vector result( numberings_[level].size(), arma::fill::zeros );
    for ( const auto& t: grids_[level] )
    {
        const math::vector elem_vec = rhs<gorder,aorder>(t);
        const auto& dofs = g_N.get_dofs( t );

        for ( size_t i = 0; i < elem_vec.n_rows; ++i )
            result( numberings_[ level ][ dofs[ i ] ] ) -= elem_vec( i );
    }

    return convert(level,result);
}

template <uint eorder, uint gorder, uint aorder>
auto solver<eorder,gorder,aorder>::convert( const fct& f ) -> math::vector
{
    unsigned char level = f.get_level();
    math::vector result( numberings_[level].size(), arma::fill::zeros );

    for ( const auto& n: numberings_[level] )
    {
        point pos = n.first;
        size_t num = n.second;
        result( num ) = f( pos );
    }

    return std::move( result );
}

template <uint eorder, uint gorder, uint aorder>
auto solver<eorder,gorder,aorder>::convert( unsigned char level, const math::vector& v ) -> fct 
{
    grid_function<real,gorder,aorder> result( grids_[level] );
    for ( auto& p: result ) p.second = 0;

    for ( const auto& n: numberings_[level] )
    {
        point  pos = n.first;
        size_t num = n.second;
        result( pos ) = v( num );
    }

    return std::move( result );
}

template <uint eorder, uint gorder, uint aorder>
void solver<eorder,gorder,aorder>::solve( fct& sol, const fct& g_N )
{
    const fct rhs = compute_rhs( g_N );
    unsigned char level = rhs.get_level();

    // Initialise defect.
    if ( level != 0 )
    {
    strategies_[ level ].apply( grids_[ level ], sol, defects_[ level ] );
    defects_[ level ] -= rhs;
    }

    real residual { 0 }; size_t counter { 0 };
    do
    {
        residual = mgm( sol, rhs );
        std::cout << "Multigrid-iteration " << ++counter
                  << ". Residual: " << residual << "." << std::endl;
    }
    while ( residual > std::pow( 2., -3*level ) && counter < 20 );
}

template <uint eorder, uint gorder, uint aorder>
real solver<eorder,gorder,aorder>::mgm( fct& sol, const fct& rhs )
{
    constexpr size_t nu1 = 3;    // Number of  pre-smoothing iterations.
    constexpr size_t nu2 = 0;    // Number of post-smoothing iterations.

    unsigned char level = rhs.get_level();

    if ( level == 0 )
    {
        cg( sol, rhs );
        return 0;
    }

    auto& S = strategies_[ level ];

    // These values will get overwritten, but copy construction
    // should be way more faster than construction from scratch.
    fct d_coarse = defects_[ level - 1 ];
    fct c_coarse = defects_[ level - 1 ];

    for ( size_t i { 0 }; i < nu1; ++i )
        smooth( sol, rhs );

    prols_[ level - 1 ].restrict( d_coarse, defects_[ level ] );
    
    // We don’t have a meaningful starting value for c_coarse,
    // starting from zero allows us to easily initialise the
    // lower level defect vector.
    defects_[ level - 1 ]  = c_coarse = 0;
    defects_[ level - 1 ] -= d_coarse;

    mgm( c_coarse, d_coarse );

    fct c_fine = defects_[ level ];
    prols_[ level - 1 ].prolongate( c_coarse, c_fine );
    sol -= c_fine;

    // Reinitialise defect.
    S.apply( grids_[ level ], sol, defects_[ level ] );
    real residual { 0 };
    for ( auto& val: defects_[ level ] )
    {
        val.second -= rhs( val.first );
        residual += val.second*val.second;
    }
    residual = std::sqrt(residual);

    for ( size_t i { 0 }; i < nu2 && residual < 1e-9; ++i )
        residual = smooth( sol, rhs );

    return residual;
}

template <uint eorder, uint gorder, uint aorder>
real solver<eorder,gorder,aorder>::smooth( fct& sol, const fct& rhs )
{
    unsigned char level = rhs.get_level();

    auto& S    = strategies_[ level ];
    auto& numb = numberings_[ level ];
    auto& diag = diagonals_ [ level ];
    auto& d    = defects_   [ level ];

    for ( auto& p: sol )
        p.second -= (2./3.) * ( d(p.first) / diag(numb[p.first])  );

    S.apply( grids_[ level ], sol, d );
    real residual { 0 };
    for ( auto& val: d ) 
    {
        val.second -= rhs( val.first );
        residual += val.second*val.second;
    }
    return std::sqrt( residual );
}

template <uint eorder, uint gorder, uint aorder>
void solver<eorder,gorder,aorder>::cg( fct& sol, const fct& rhs )
{
    const unsigned char lvl = rhs.get_level();
    auto& S    = strategies_[ lvl ];
    auto& diag =  diagonals_[ lvl ];

    const math::vector b = convert( rhs );
    math::vector x = convert( sol );


    fct zz(sol);
    fct dd(sol);

    real alpha, beta;
    math::vector r, h, d, z;

    r.zeros( b.size() );
    h.zeros( b.size() );
    d.zeros( b.size() );
    z.zeros( b.size() );

    S.apply( grids_[ lvl ], dd, zz );
    for ( auto& val: zz )
    {
        val.second = rhs( val.first ) - val.second;
    }
    
    r = convert( zz );
    h = r / diag;

     d = h;
    dd = convert( lvl, d );

    for ( size_t k = 0; k <= 100; ++k )
    {
        S.apply( grids_[ lvl ], dd, zz );
        z = convert( zz );

        alpha = dot( r, h ) / dot( d, z );

        math::vector rold = r;
        math::vector hold = h;

        x = x + alpha*d;
        r = r - alpha*z;
        h = r / diag;

        if ( norm(r) < 1e-9 || k == 100 )
            break;

        beta = dot( r, h )/dot( rold, hold );
        d = h + beta*d;
        dd = convert( lvl, d );
    }
    sol = convert( lvl, x );
}

}

