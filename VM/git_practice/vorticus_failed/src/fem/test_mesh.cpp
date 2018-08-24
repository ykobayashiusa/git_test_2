/*
 * Copyright (C) 2016 Matthias Kirchhart
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

#include "fem/multigrid.h"
#include "fem/vtkwriter.h"
#include "fem/gmsh_reader.h"

#include <fstream>
#include <sstream>
#include <iostream>

using namespace fem;

int main( int argc, char *argv[] )
{
    using std::cout;
    if ( argc < 2 )
    {
        std::cerr << "Expected filename as parameter.\n";
        return -1;
    }
    
    multigrid<1> mg = read_gmsh( argv[1] );

    for ( size_t lvl = 0; lvl < 6; ++lvl )
    {
        cout << "Refining mesh..."; cout.flush();
        mg.adapt();
        cout << "done.\n"; cout.flush();


        cout << "Numbering mesh..."; cout.flush();
        node_numbering<1,1>      numb( mg.grid_begin(lvl), mg.grid_end(lvl) );
        cout << "done.\n"; cout.flush();
        grid_function<real,1,1>  numb_fct( numb );

        for ( size_t i = 0; i < numb.size(); ++i )
            numb_fct( numb(i) ) = i;


        vtkwriter<1> writer( mg.grid_begin(lvl), mg.grid_end(lvl) );
        writer.register_scalar( numb_fct, "number" );

        std::stringstream filename; filename << "lol_" << lvl << ".vtu";
        std::ofstream file( filename.str() );
        writer.write( file );
        cout << "Marking mesh..."; cout.flush();
        for ( auto it = mg.grid_begin(lvl); it != mg.grid_end(lvl); ++it )
            it->set_ref();
        cout << "done.\n"; cout.flush();
    }
}

