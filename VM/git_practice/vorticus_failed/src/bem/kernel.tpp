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

#include "bem/differentials.h"

namespace bem
{

/*!
 * \brief Evaluates the regularised hypersingular kernel W.
 * \param t1  The ‘outer’ x-triangle
 * \param pos The position on t1.
 * \param t2  The ‘inner’ y-triangle
 * \param pos The position on t2.
 *
 * This function evaluates the regularised kernel function for the
 * hypersingular kernel W and all pairs of shapefunctions of the
 * involved triangles.
 * \f[
 * \frac{\langle \mbox{rot}_\Gamma\, N_1(\mathbf{x}),\:
 *               \mbox{rot}_\Gamma\, N_2(\mathbf{y})\rangle}
 *      {4\pi\vert\mathbf{x}-\mathbf{y}\vert}
 * \f]
 */
template <uint gorder, uint aorder> inline
auto W_kernel( const triangle<gorder>& t1, const bary pos1,
               const triangle<gorder>& t2, const bary pos2 )
-> arma::mat::fixed< shapefcts2d::num<aorder>(), shapefcts2d::num<aorder>() >
{
    using math::pifac;
    using geometry::point;

    const point  p1 = t1.Chi      ( pos1 );
    const auto rot1 = rotN<gorder,aorder>( t1, pos1 );

    const point  p2 = t2.Chi      ( pos2 );
    const auto rot2 = rotN<gorder,aorder>( t2, pos2 );

    const real r = ( p1 - p2 ).r();
    constexpr size_t size = shapefcts2d::num<aorder>();
    arma::mat::fixed<size,size> result;
    for ( size_t i { 0 }; i < size; ++i )
    {
        for ( size_t j { 0 }; j < size; ++j )
        {
            result(i,j) = (pifac/r) * scal_prod( rot1[i], rot2[j] );
        }
    }

    return result;
}

}

