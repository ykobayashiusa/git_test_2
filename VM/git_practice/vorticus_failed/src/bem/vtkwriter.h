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
#ifndef BEM_VTK_WRITER_H
#define BEM_VTK_WRITER_H

#include <map>
#include <string>
#include <ostream>
#include <unordered_map>

#include "bem/multigrid.h"
#include "bem/grid_function.h"

namespace bem
{

template <uint order>
class vtkwriter
{
public:
    vtkwriter( const typename multigrid<order>::grid& triang );

    void write( std::ostream& os );

    using string = std::string;
    using scalar_function = grid_function_base<real, order>;
    using vector_function = grid_function_base<point,order>;
    void register_scalar( string name, const scalar_function& f );
    void register_vector( string name, const vector_function& f );

private:
    const typename multigrid<order>::grid&  triang_;
    std::unordered_map<point,size_t>                 node_numbers_;
    std::unordered_map<size_t,point>                 node_numbers_inv_;

    std::map<string,const scalar_function*> scalars_;
    std::map<string,const vector_function*> vectors_;
};

}

#include "bem/vtkwriter.tpp"
#endif

