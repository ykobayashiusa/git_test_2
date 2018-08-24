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
#ifndef VRM_PANEL_STRATEGY_H
#define VRM_PANEL_STRATEGY_H

#include "fmm/fmm.h"
#include "math/constants.h"
#include "math/harmonics.h"

namespace vrm
{

template <uint eorder, typename accessor>
class particle_strategy
{
public:
    using point       = geometry::point;
    using cmplx_point = geometry::cmplx_point;
    using Mcoeff_t    = std::array<cmplx,eorder>;
    using Lcoeff_t    = std::array<cmplx,eorder>;

    particle_strategy();
    particle_strategy( accessor p_get ): get(std::move(p_get)) {}

    template <typename source_iterator>
    void p2m( Mcoeff_t &M, point xc, source_iterator *const begin,
                                     source_iterator *const end ) const noexcept;
    template <typename source_iterator>
    void p2l( Lcoeff_t &L, point xc, source_iterator *const begin,
                                     source_iterator *const end ) const noexcept;

    void m2m( Mcoeff_t &target, const Mcoeff_t &source, point r ) const noexcept;
    void m2l( Lcoeff_t &target, const Mcoeff_t &source, point r ) const noexcept;
    void l2l( Lcoeff_t &target, const Lcoeff_t &source, point r ) const noexcept;

    template <typename target_iterator>
    void m2p( target_iterator *const begin, target_iterator *const end,
              const Mcoeff_t& M, point xc ) noexcept;

    template <typename target_iterator>
    void l2p( target_iterator *const begin, target_iterator *const end,
              const Lcoeff_t& M, point xc ) noexcept;

    template <typename target_iterator, typename source_iterator>
    void p2p( target_iterator *const tbegin, target_iterator *const tend,
              source_iterator *const sbegin, source_iterator *const send ) noexcept;

    bool mac( point Acentre, real Asize, point Bcentre, real Bsize ) const noexcept;

    uint source_leaf_max() const noexcept { return 16; }
    uint target_leaf_max() const noexcept { return 16; }

    static constexpr real time_m2l()             noexcept { return 0; }
    static constexpr real time_m2p( uint )       noexcept { return 1; }
    static constexpr real time_p2l( uint )       noexcept { return 1; }
    static constexpr real time_p2p( uint, uint ) noexcept { return 1; }

    template <typename iterator>
    geometry::point pos( iterator it ) const noexcept;

    template <typename iterator>
    void   bounding_box( iterator it, geometry::point& min, geometry::point& max ) const noexcept;

    std::array<cmplx,eorder> I( point x ) const noexcept;
    std::array<cmplx,eorder> O( point x ) const noexcept;

    std::array<cmplx,eorder> dI( point x ) const noexcept;
    std::array<cmplx,eorder> dO( point x ) const noexcept;

    accessor get {};
    real    epsilon;
    real    theta_max { 0.8 };
};

}

#include "vrm/particle_strategy.tpp"
#endif

