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
 * \brief Computes the surface-curl of the shape-functions on a triangle.
 * \param t   The triangle on which the shape-functions are defined.
 * \param pos The position at which the functions should be evaluated.
 *
 * This function computes the curl of the shape-functions of order aorder
 * on the given triangle. It is defined by
 * \f[
 * \mbox{rot}(N(\mathbf{x})\mathbf{n}(\mathbf{x})) = \mbox{rot}_\Gamma\bigl(N(\mathbf{x})\bigr) =
 * \mathbf{J}(\mathbf{J}^{\top}\mathbf{J} )^{-1} N(\mathbf{x}) \times \mathbf{n},
 * \f]
 *
 * where \f$\mathbf{J}\f$ denotes the Jacobian of the mapping \f$\chi\f$ from
 * the reference triangle to space. As it depends on the particular triangle,
 * it cannot be defined in the file shape_functions.h.
 */
template <uint gorder, uint aorder> inline
shapefcts2d::coeffs<aorder,point> rotN( const triangle<gorder> &t, const bary pos ) noexcept
{
    const point  xi  = t.dChidxi ( pos );
    const point  eta = t.dChideta( pos );

    const auto N_xi  = shapefcts2d::dNdxi <aorder>( pos );
    const auto N_eta = shapefcts2d::dNdeta<aorder>( pos );

    const point un      = cross_prod( xi, eta );
    const real  detG    =  scal_prod( un, un  );
    const point n       = un/std::sqrt(detG);

    const real    J[3][2] = { { xi.x, eta.x },
                              { xi.y, eta.y },
                              { xi.z, eta.z } };
    const real invG[2][2] = { {  scal_prod(eta,eta)/detG, -scal_prod( xi, eta)/detG },
                              { -scal_prod( xi,eta)/detG,  scal_prod( xi,  xi)/detG } };

    const real map[3][2]  = { { J[0][0]*invG[0][0] + J[0][1]*invG[1][0],    J[0][0]*invG[0][1] + J[0][1]*invG[1][1] },
                              { J[1][0]*invG[0][0] + J[1][1]*invG[1][0],    J[1][0]*invG[0][1] + J[1][1]*invG[1][1] },
                              { J[2][0]*invG[0][0] + J[2][1]*invG[1][0],    J[2][0]*invG[0][1] + J[2][1]*invG[1][1] } };


    shapefcts2d::coeffs<aorder,point> result;
    for ( size_t i { 0 }; i < result.size(); ++i )
    {
        result[ i ] = cross_prod( point { map[0][0]*N_xi[i] + map[0][1]*N_eta[i],
                                          map[1][0]*N_xi[i] + map[1][1]*N_eta[i],
                                          map[2][0]*N_xi[i] + map[2][1]*N_eta[i] },
                                  n );
    }

    return result;
}

/*!
 * \brief Computes the curl of the shape-functions on a triangle.
 * \param t   The triangle on which the shape-functions are defined.
 * \param pos The position at which the functions should be evaluated.
 */
template <uint gorder, uint aorder> inline
auto rotNvec( const triangle<gorder> &t, const bary pos ) noexcept
-> shapefcts2d::coeffs<aorder,std::array<point,3>>
{
    const point  xi  = t.dChidxi ( pos );
    const point  eta = t.dChideta( pos );

    const auto N_xi  = shapefcts2d::dNdxi <aorder>( pos );
    const auto N_eta = shapefcts2d::dNdeta<aorder>( pos );

    const point un      = cross_prod( xi, eta );
    const real  detG    =  scal_prod( un, un  );

    const real    J[3][2] = { { xi.x, eta.x },
                              { xi.y, eta.y },
                              { xi.z, eta.z } };
    const real invG[2][2] = { {  scal_prod(eta,eta)/detG, -scal_prod( xi, eta)/detG },
                              { -scal_prod( xi,eta)/detG,  scal_prod( xi,  xi)/detG } };

    const real map[3][2]  = { { J[0][0]*invG[0][0] + J[0][1]*invG[1][0],    J[0][0]*invG[0][1] + J[0][1]*invG[1][1] },
                              { J[1][0]*invG[0][0] + J[1][1]*invG[1][0],    J[1][0]*invG[0][1] + J[1][1]*invG[1][1] },
                              { J[2][0]*invG[0][0] + J[2][1]*invG[1][0],    J[2][0]*invG[0][1] + J[2][1]*invG[1][1] } };


    shapefcts2d::coeffs<aorder,std::array<point,3>> result;
    for ( size_t i { 0 }; i < result.size(); ++i )
    {
        const point grad =  point { map[0][0]*N_xi[i] + map[0][1]*N_eta[i],
                                    map[1][0]*N_xi[i] + map[1][1]*N_eta[i],
                                    map[2][0]*N_xi[i] + map[2][1]*N_eta[i] };


        result[ i ][ 0 ] = point { 0, grad.z, -grad.y };
        result[ i ][ 1 ] = point { -grad.z, 0, grad.x };
        result[ i ][ 2 ] = point { grad.y, -grad.x, 0 };
    }

    return result;
}

}

