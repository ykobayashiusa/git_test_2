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

#include "bem/dunavant.h"
#include "bem/quadrature.h"
#include "bem/differentials.h"

namespace bem
{

template <uint gorder, uint aorder>
auto sheet_functional( const grid<gorder>& g,
                       const grid_function<real,gorder,aorder>& mu )
-> grid_function<point,gorder,aorder>
{
    grid_function<point,gorder,aorder> result(g);

    constexpr uint  quad_order   = 2*(aorder + 1);
    const     uint  quad_num     = dunavant_num_points[ quad_order ];
    const     bary* quad_points  = dunavant_points    [ quad_order ];
    const     real* quad_weights = dunavant_weights   [ quad_order ];

    for ( const triangle<gorder>& t: g )
    {
        shapefcts2d::coeffs<aorder,point> triangle_results;
        for ( uint q = 0; q < quad_num; ++q )
        {
            const bary pos = quad_points [ q ];
            const real w   = quad_weights[ q ];
            const real gx  = t.surf_elem( pos );
          
            const auto  rot  = rotNvec<gorder,aorder>( t, pos );
            const point mu_n = t.normal( pos ) * mu.eval( t, pos );

            for ( size_t j = 0; j < triangle_results.size(); ++j )
            {
                triangle_results[ j ] += -0.5*gx*w*point { scal_prod( mu_n, rot[ j ][ 0 ] ),
                                                           scal_prod( mu_n, rot[ j ][ 1 ] ),
                                                           scal_prod( mu_n, rot[ j ][ 2 ] ) };
            }
        }

        const auto& dof_map = result.get_dofs( t );
        for ( size_t j = 0; j < triangle_results.size(); ++j )
        {
            result( dof_map[ j ] ) += triangle_results[ j ];
        }
    }
    return std::move(result);
}

template <uint gorder, uint aorder>
grid_function<point,gorder,aorder> compute_sheet( const grid<gorder>& g,
                                                  const grid_function<real,gorder,aorder>& mu )
{
    const grid_function<point,gorder,aorder> rhs    = sheet_functional( g, mu );
          grid_function<point,gorder,aorder> result( g ); result = point {};

    const sparse_matrix<gorder,aorder> mass = mass_matrix<gorder,aorder>( g );

    grid_function<point,gorder,aorder> r = mass.apply( result );
    r -= rhs; r *= -1;
    
    grid_function<point,gorder,aorder> d = r;
    
    real residual = 0;
    size_t counter = 0;
    do
    {
        auto z = mass.apply( d );

        real tmp1 = 0;
        for ( const auto& val: r ) tmp1 += scal_prod( val.second, val.second );
        real tmp2 = 0;
        for ( const auto& val: d ) tmp2 += scal_prod( val.second, z( val.first ) );
        real alpha = tmp1/tmp2;

        auto tmp3 = d; tmp3*=alpha; result += tmp3;
        tmp3 = z; tmp3 *= alpha; r -= tmp3;

        residual = 0;
        for ( const auto& val: r ) residual += scal_prod( val.second, val.second );
        residual = std::sqrt( residual );
    
        std::cout << "Iteration " << counter++ << ". Residual " << residual << ".\n";
        
        tmp2 = 0;
        for ( const auto& val: r ) tmp2 += scal_prod( val.second, val.second );
        real beta = tmp2 / tmp1;

        d *= beta; d+= r;
    }
    while ( residual > 1e-9 && counter < 10000 );

    return std::move( result );
}

}

