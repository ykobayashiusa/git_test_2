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

#include "bem/solver.h"
#include "bem/vtkwriter.h"
#include "misc/stopwatch.h"
#include "bem/sheet.h"

#include <fstream>
#include <sstream>
#include <iostream>


using namespace bem;

constexpr unsigned int  gorder = 1;
constexpr unsigned int  aorder = 1;
constexpr unsigned int  maxlvl = 6;

void refine( typename multigrid<gorder>::grid &g );
void  solve( multigrid<gorder>& mg );

int main()
{
    multigrid<gorder> mg;

    for ( unsigned int l = 0; l < maxlvl; ++l )
    {
        typename multigrid<gorder>::grid g(mg,l);
        refine(g);
    }

    solve( mg );
}

void refine( typename multigrid<gorder>::grid &g )
{
    for ( triangle<gorder>& t: g )
    {
        t.set_ref_mark();
    }
    g.get_multi_grid().adapt();
}

void solve( multigrid<gorder>& mg )
{
    std::vector<typename multigrid<gorder>::grid> grids;
    std::vector<grid_function<double,gorder,aorder>> fcts;
    std::vector<prolongation<gorder,aorder>> prols;

    for ( size_t l = 0; l <= mg.last_level(); ++l )
    {
        grids.push_back(typename multigrid<gorder>::grid(mg,l));
        auto& g = grids.back();
        fcts.push_back(grid_function<double,gorder,aorder>( g ) );
        if ( l >= 1 ) prols.push_back( prolongation<gorder,aorder>( grids[ l - 1 ], grids[ l ] ) );
        std::cout << "Level " << l << ": triangles: "
                  << g.size() << ", dofs: "
                  << fcts.back().size() << "." << std::endl;
    }

    std::cout << "Creating solver…"; std::cout.flush();
    solver<6,gorder,aorder> solv(mg);
    std::cout << "done.\n"; std::cout.flush();

    for ( size_t l = 0; l <= mg.last_level(); ++l )
    {
        std::cout << "Solving on level " << l << std::endl;
        typename multigrid<gorder>::grid g(mg,l);
        
        if ( l >= 1 ) prols[ l - 1 ].prolongate( fcts[ l - 1 ], fcts[ l ] );

        misc::stopwatch watch;
        solv.solve( fcts[ l ], fcts[ l ] );
        std::cout << "done after " << watch.elapsed()
                  << " seconds."   << std::endl;

        auto test = compute_sheet( grids[ l ], fcts[ l ] );

        vtkwriter<gorder> w( g );
        w.register_scalar( "mu", fcts[ l ] );
        w.register_vector( "omega", test );

        std::stringstream filename;
        filename << "mu_" << l << ".vtu";
        std::ofstream file( filename.str() );
        w.write( file );
    }
    for ( size_t l = 0; l <= mg.last_level(); ++l )
    {
        std::cout << "Solving on level " << l << std::endl;
        typename multigrid<gorder>::grid g(mg,l);
        
        if ( l >= 1 ) prols[ l - 1 ].prolongate( fcts[ l - 1 ], fcts[ l ] );

        misc::stopwatch watch;
        solv.solve( fcts[ l ], fcts[ l ] );
        std::cout << "done after " << watch.elapsed()
                  << " seconds."   << std::endl;

        auto test = compute_sheet( grids[ l ], fcts[ l ] );

        vtkwriter<gorder> w( g );
        w.register_scalar( "mu", fcts[ l ] );
        w.register_vector( "omega", test );

        std::stringstream filename;
        filename << "mu_" << l << ".vtu";
        std::ofstream file( filename.str() );
        w.write( file );
    }
}

