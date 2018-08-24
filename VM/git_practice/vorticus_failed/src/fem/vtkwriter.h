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
#ifndef FEM_VTK_WRITER_H
#define FEM_VTK_WRITER_H

#include <map>
#include <string>
#include <ostream>
#include <unordered_map>

#include "fem/multigrid.h"
#include "fem/grid_function.h"

namespace fem
{

template <size_t order>
class vtkwriter
{
public:
    vtkwriter( const_grid_iterator<order> begin, const_grid_iterator<order> end );
    void write( std::ostream& os );

    using scalar = grid_function_base<real,order>;
    void register_scalar( const scalar &f, std::string name );

    using vector = grid_function_base<point,order>;
    void register_vector( const vector &f, std::string name );

private:
    node_numbering<order,1> numb_;

    std::map<std::string,const scalar*> scalars_;
    std::map<std::string,const vector*> vectors_;
};

}

#include "fem/vtkwriter.tpp"
#endif

