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
#ifndef FEM_GMSH_READER_H
#define FEM_GMSH_READER_H

#include <string>
#include "fem/multigrid.h"
#include "fem/multigrid_builder.h"

namespace fem
{

class gmsh_reader: public multigrid_builder<1>
{
public:
    gmsh_reader( const std::string p_filename ): filename { p_filename } {}
    multigrid<1> make_grid();

private:
    std::string filename;
};

inline
multigrid<1> read_gmsh( const std::string &filename )
{
    return gmsh_reader( filename ).make_grid();
}

}

#endif

