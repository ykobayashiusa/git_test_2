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
#ifndef BEM_PANEL_STRATEGY_H
#define BEM_PANEL_STRATEGY_H

#include "bem/kernel.h"
#include "bem/dunavant.h"
#include "bem/multigrid.h"
#include "bem/quadrature.h"
#include "bem/grid_function.h"

#include "fmm/fmm.h"
#include "math/harmonics.h"

namespace bem
{

template <uint eorder, uint gorder, uint aorder>
class panel_strategy
{
public:
    static constexpr size_t tdofs { lattice::size<aorder>() };

    using point       = geometry::point;
    using cmplx_point = geometry::cmplx_point;
    using grid_t      = typename multigrid<gorder>::grid;
    using source_iter = typename grid_t::const_iterator;
    using target_iter = typename grid_t::const_iterator;
    using Mcoeff_t    = math::cmplx_point_coeff<eorder>;
    using Lcoeff_t    = math::cmplx_point_coeff<eorder>;

    using elem_mat_t = arma::mat::fixed< shapefcts2d::num<aorder>(), shapefcts2d::num<aorder>() >;

    panel_strategy( const grid_t& g );
                   

    void apply( const grid_t& grid,
                const grid_function<real,gorder,aorder>& source, grid_function<real,gorder,aorder>& target );

    void p2m( Mcoeff_t &M, point xc, const source_iter *const begin,
                                     const source_iter *const end ) const noexcept;
    void p2l( Mcoeff_t &M, point xc, const source_iter *const begin,
                                     const source_iter *const end ) const noexcept;

    void m2m( Mcoeff_t &target, const Mcoeff_t &source, point r ) const noexcept;
    void m2l( Lcoeff_t &target, const Mcoeff_t &source, point r ) const noexcept;
    void l2l( Lcoeff_t &target, const Lcoeff_t &source, point r ) const noexcept;

    void m2p( target_iter *const begin, target_iter *const end,
              const Mcoeff_t& M, point xc ) noexcept;

    void l2p( target_iter *const begin, target_iter *const end,
              const Lcoeff_t& M, point xc ) noexcept;

    void p2p(       target_iter *const tbegin,       target_iter *const tend,
              const source_iter *const sbegin, const source_iter *const send ) noexcept;

    bool mac( point Acentre, real Asize, point Bcentre, real Bsize ) const noexcept;

    uint source_leaf_max() const noexcept { return 1; }
    uint target_leaf_max() const noexcept { return 1; }

    const elem_mat_t& get_elem_mat( const triangle<gorder>& t1, const triangle<gorder>& t2 );

    static constexpr real time_m2l()             noexcept { return 0; }
    static constexpr real time_p2l( uint )       noexcept { return 1; }
    static constexpr real time_m2p( uint )       noexcept { return 1; }
    static constexpr real time_p2p( uint, uint ) noexcept { return 1; }

    static geometry::point pos( source_iter it ) noexcept
    {
        return it->pos();
    }

    static void bounding_box( source_iter it, geometry::point& min, geometry::point& max ) noexcept
    {
        it->bounding_box( min, max );
    }

private:
    point eval_l_exp( const Lcoeff_t &L, point r ) const noexcept;
    point eval_m_exp( const Lcoeff_t &L, point r ) const noexcept;

    std::unordered_map<point, std::unordered_map<point,elem_mat_t> > elem_matrices;
    std::unordered_map<geoid, shapefcts2d::vals<aorder> >            partial_results;

    const grid_function<real,gorder,aorder>* source_strengths { nullptr };
};

}

#include "panel_strategy.tpp"
#endif

