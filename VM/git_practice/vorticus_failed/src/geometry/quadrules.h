/*
 * Copyright (C) 2017 Matthias Kirchhart
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
#ifndef GEOMETRY_QUADRULES_H
#define GEOMETRY_QUADRULES_H

#include "geometry/bary.h"
#include "geometry/point.h"

#include <vector>
#include <utility>

namespace geometry
{

struct quad_node
{
    point x;
    real  w;
};

struct triangular_quad_node
{
    bary2d b;
    real   w;
};

struct tetrahedral_quad_node
{
    bary3d b;
    real   w;
};

using  triangular_quadrule = std::vector< triangular_quad_node>;
using tetrahedral_quadrule = std::vector<tetrahedral_quad_node>;
using     cubical_quadrule = std::vector<quad_node>;

const  triangular_quadrule&  get_triangular_quadrule( size_t degree );
const tetrahedral_quadrule& get_tetrahedral_quadrule( size_t degree );
const     cubical_quadrule      get_cubical_quadrule( size_t degree );

tetrahedral_quadrule get_refined_tetrahedral_quadrule( size_t degree, size_t refinements );
    cubical_quadrule get_refined_cubical_quadrule    ( size_t degree, size_t refinements );

}

#endif

