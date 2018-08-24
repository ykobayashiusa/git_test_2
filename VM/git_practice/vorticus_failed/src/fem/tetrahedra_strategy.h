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
#ifndef FEM_TETRAHEDRA_STRATEGY_H
#define FEM_TETRAHEDRA_STRATEGY_H

#include "fem/multigrid.h"
#include "fem/potentials.h"
#include "fem/grid_function.h"

#include "fmm/fmm.h"
#include "math/harmonics.h"
#include "math/constants.h"
#include "geometry/quadrules.h"

namespace fem
{

template <uint eorder>
class tetrahedra_strategy
{
public:
    using point       = geometry::point;
    using cmplx_point = geometry::cmplx_point;
    using Mcoeff_t    = math::cmplx_point_coeff<eorder>;
    using Lcoeff_t    = math::cmplx_point_coeff<eorder>;

    const grid_function<point,1,1>* omega { nullptr };
    
    template <typename source_iter>               
    void p2m( Mcoeff_t &M, point xc, const source_iter *const begin,
                                     const source_iter *const end ) const;
    template <typename source_iter>               
    void p2l( Mcoeff_t &M, point xc, const source_iter *const begin,
                                     const source_iter *const end ) const;

    void m2m( Mcoeff_t &target, const Mcoeff_t &source, point r ) const noexcept;
    void m2l( Lcoeff_t &target, const Mcoeff_t &source, point r ) const noexcept;
    void l2l( Lcoeff_t &target, const Lcoeff_t &source, point r ) const noexcept;

    template <typename target_iter>
    void m2p( target_iter *const begin, target_iter *const end,
              const Mcoeff_t& M, point xc ) const;

    template <typename target_iter>
    void l2p( target_iter *const begin, target_iter *const end,
              const Lcoeff_t& M, point xc ) const;

    template <typename target_iter, typename source_iter>
    void p2p(       target_iter *const tbegin,       target_iter *const tend,
              const source_iter *const sbegin, const source_iter *const send ) const;

    bool mac( point Acentre, real Asize, point Bcentre, real Bsize ) const noexcept;

    uint source_leaf_max() const noexcept { return  8; }
    uint target_leaf_max() const noexcept { return 32; }

    real time_m2l()             const noexcept { return 0;    }
    real time_p2l( uint )       const noexcept { return 1;    }
    real time_m2p( uint )       const noexcept { return 1;    }
    real time_p2p( uint, uint ) const noexcept { return 1e20; }

    template <typename source_iter>
    geometry::point pos( source_iter it ) const noexcept;

    template <typename source_iter>
    void bounding_box( source_iter it, geometry::point& min, geometry::point& max ) const noexcept;
};

}

#include "fem/tetrahedra_strategy.tpp"
#endif

