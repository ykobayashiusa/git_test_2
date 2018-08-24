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
#ifndef FEM_NODE_NUMBERING_H
#define FEM_NODE_NUMBERING_H

#include "fem/multigrid.h"
#include "geometry/zsort.h"

namespace fem
{

template <size_t gorder, size_t aorder>
class node_numbering
{
public:
    using coords_t = shapefcts3d::coeffs<aorder,point>;
    using  numbs_t = shapefcts3d::coeffs<aorder,size_t>;

    node_numbering( const_grid_iterator<gorder> begin, const_grid_iterator<gorder> end );

    node_numbering() = default;
    node_numbering( const node_numbering&  ) = default;
    node_numbering(       node_numbering&& ) = default;
    node_numbering& operator=( const node_numbering&  ) = default;
    node_numbering& operator=(       node_numbering&& ) = default;
    ~node_numbering() = default;

    const std::vector<point>&                      nodes() const noexcept { return      nodes_; }
    const std::unordered_map<point,size_t>&    inv_nodes() const noexcept { return  inv_nodes_; }
    const std::unordered_map<geoid,numbs_t>&  cell_nodes() const noexcept { return cell_nodes_; }

    size_t size() const noexcept { return nodes_.size(); }

    size_t  operator()( point  x ) const;
    point   operator()( size_t n ) const;
    numbs_t operator()( const tetrahedron<gorder> &t ) const;

private:
    static coords_t get_dofs( const tetrahedron<gorder> &t ) noexcept;
    static point    get_face_dof( const tetrahedron<gorder> &t, bary3d b ) noexcept;
    static point    get_edge_dof( const tetrahedron<gorder> &t, bary3d b ) noexcept;
    static point    get_node_dof( const tetrahedron<gorder> &t, bary3d b ) noexcept;

private:
    std::vector<point>                     nodes_;
    std::unordered_map<point,size_t>   inv_nodes_;
    std::unordered_map<geoid,numbs_t> cell_nodes_;
};

}

#include "fem/node_numbering.tpp"
#endif

