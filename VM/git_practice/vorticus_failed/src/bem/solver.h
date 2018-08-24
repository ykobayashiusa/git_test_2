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
#ifndef BEM_SOLVER_H
#define BEM_SOLVER_H

#include "bem/multigrid.h"
#include "bem/panel_strategy.h"

#include <vector>
#include <memory>
#include <unordered_map>

namespace bem
{

template <uint eorder, uint gorder, uint aorder>
class solver
{
public:
    solver( multigrid<gorder>& mg );

    using fct = grid_function<real,gorder,aorder>;
    void solve( fct& sol, const fct& g_N );
    
private:
    fct          convert( unsigned char level, const math::vector& x );
    math::vector convert( const fct& f );

    fct compute_rhs( const fct& g_N );

    void     cg( fct& sol, const fct& rhs );
    real    mgm( fct& sol, const fct& rhs );
    real smooth( fct& sol, const fct& rhs );

private:
    multigrid<gorder>& mg_;
    using grid     = typename multigrid<gorder>::grid;
    using strategy = panel_strategy<eorder,gorder,aorder>;
    using prol     = prolongation<gorder,aorder>;
    
    std::vector<grid>      grids_;
    std::vector<strategy>  strategies_;

    std::vector< std::unordered_map<point,size_t> >     numberings_;
    std::vector< math::vector >                         diagonals_;

    std::vector< fct >  defects_;
    std::vector< prol > prols_;
};

}

#include "bem/solver.tpp"
#endif

