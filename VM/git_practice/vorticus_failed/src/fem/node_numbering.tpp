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

namespace fem
{

template <size_t gorder,size_t aorder>
node_numbering<gorder,aorder>::node_numbering
(
    const_grid_iterator<gorder> begin,
    const_grid_iterator<gorder> end
)
{
    // Generate all nodes.
    for ( auto it = begin; it != end; ++it )
    {
        const coords_t dofs { get_dofs(*it) };
        for ( point dof: dofs ) inv_nodes_[ dof ];
    }
 
    // Sort nodes according to Z-order. 
    nodes_.reserve( inv_nodes_.size() );
    for ( auto it = inv_nodes_.begin(); it != inv_nodes_.end(); ++it )
    {
        nodes_.push_back( it->first );
    }
    geometry::zsort( nodes_.begin(), nodes_.end() );
   
    // Generate inverse of the numbering. 
    for ( size_t i = 0; i < nodes_.size(); ++i )
    {
        inv_nodes_[ nodes_[ i ] ] = i;
    }
    
    // Store numbers of nodes associated with cells.
    for ( auto it = begin; it != end; ++it )
    {
        const coords_t dofs { get_dofs(*it) };
        numbs_t result;
        for ( size_t i = 0; i < dofs.size(); ++i )
        {
            result[i] = inv_nodes_[dofs[i]];
        }
        cell_nodes_[ it->id() ] = result;
    }
}

template <size_t gorder, size_t aorder> inline
size_t node_numbering<gorder,aorder>::operator()( point x ) const
{
    auto it = inv_nodes_.find(x);
    if ( it == inv_nodes_.end() )
        throw std::out_of_range { "Node numbering: point -> size_t." };
    return it->second;
}

template <size_t gorder, size_t aorder> inline
point node_numbering<gorder,aorder>::operator()( size_t n ) const
{
    if ( n >= nodes_.size() )
        throw std::out_of_range { "Node numbering: size_t -> point." };
    return nodes_[ n ];
}

template <size_t gorder, size_t aorder> inline
auto node_numbering<gorder,aorder>::operator()( const tetrahedron<gorder> &t ) const 
-> numbs_t
{
    auto it = cell_nodes_.find(t.id());
    if ( it == cell_nodes_.end() )
        throw std::out_of_range { "Node numbering: tet -> numbs_t." };
    return it->second;
}

template <size_t gorder, size_t aorder>
auto node_numbering<gorder,aorder>::get_dofs( const tetrahedron<gorder> &t ) noexcept
-> coords_t
{
    constexpr shapefcts3d::coeffs<aorder,bary3d> lat { shapefcts3d::lattice<aorder>() };

    shapefcts3d::coeffs<aorder,point> result;
    for ( size_t i = 0; i < lat.size(); ++i )
    {
        bary3d b = lat[i];
        size_t zero_count = 0;
        if ( b.z0 == 0 ) ++zero_count;
        if ( b.z1 == 0 ) ++zero_count;
        if ( b.z2 == 0 ) ++zero_count;
        if ( b.z3 == 0 ) ++zero_count;

        switch ( zero_count )
        {
        case 0: result[i] = t.Chi( b ); break;
        case 1: result[i] = get_face_dof( t, b ); break;
        case 2: result[i] = get_edge_dof( t, b ); break;
        case 3: result[i] = get_node_dof( t, b ); break;
        }
    }

    return result;
}

template <size_t gorder, size_t aorder>
auto node_numbering<gorder,aorder>::get_face_dof( const tetrahedron<gorder> &t, bary3d b ) noexcept
-> point
{
    const triangle<gorder> *f;
    bary2d p;
    if ( b.z0 == 0 )
    {
        f = &(t.get_face(0));
        p.z0 = b.z1;
        p.z1 = b.z2;
        p.z2 = b.z3;
    }
    else if ( b.z1 == 0 )
    {
        f = &(t.get_face(1));
        p.z0 = b.z0;
        p.z1 = b.z2;
        p.z2 = b.z3;
    }
    else if ( b.z2 == 0 )
    {
        f = &(t.get_face(2));
        p.z0 = b.z0; p.z1 = b.z1;
        p.z2 = b.z3;
    }
    else // ( b.z3 == 0 )
    {
        f = &(t.get_face(3));
        p.z0 = b.z0;
        p.z1 = b.z1;
        p.z2 = b.z2;
    }

    return f->Chi( p );
}

template <size_t gorder, size_t aorder>
auto node_numbering<gorder,aorder>::get_edge_dof( const tetrahedron<gorder> &t, bary3d b ) noexcept
-> point
{
    const edge<gorder> *e { nullptr };
    real p;
    if ( b.z0 != 0 && b.z1 != 0 )
    {
        e = &(t.get_edge(0));
        p = b.z0;
    }
    else if ( b.z0 != 0 && b.z2 != 0 )
    {
        e = &(t.get_edge(1));
        p = b.z0;
    }
    else if ( b.z1 != 0 && b.z2 != 0 )
    {
        e = &(t.get_edge(2));
        p = b.z1;
    }
    else if ( b.z0 != 0 && b.z3 != 0 )
    {
        e = &(t.get_edge(3));
        p = b.z0;
    }
    else if ( b.z1 != 0 && b.z3 != 0 )
    {
        e = &(t.get_edge(4));
        p = b.z1;
    }
    else // ( b.z2 != 0 && b.z3 != 0 )
    {
        e = &(t.get_edge(5));
        p = b.z2;
    }
    return e->Chi(p);
}

template <size_t gorder, size_t aorder>
auto node_numbering<gorder,aorder>::get_node_dof( const tetrahedron<gorder> &t, bary3d b ) noexcept
-> point
{
         if ( b.z0 != 0 ) return t.get_node(0);
    else if ( b.z1 != 0 ) return t.get_node(1);
    else if ( b.z2 != 0 ) return t.get_node(2);
    else                  return t.get_node(3);
}

}

