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

template <size_t order>
class multigrid<order>::const_grid_iterator
{
public:
    using iter = typename std::unordered_map<geoid,tetrahedron<order>>::const_iterator;
    using size_type         = size_t;
    using difference_type   = size_t;
    using value_type        = const tetrahedron<order>;
    using reference         = const tetrahedron<order>&;
    using pointer           = const tetrahedron<order>*;
    using iterator_category = std::forward_iterator_tag;


    iter  it_;
    uchar level_;

    bool operator==( const const_grid_iterator &rhs ) const noexcept
    { return level_ == rhs.level_ && it_ == rhs.it_; }

    bool operator!=( const const_grid_iterator &rhs ) const noexcept
    { return level_ != rhs.level_ || it_ != rhs.it_; }

    reference operator* () const noexcept { return   it_->second; }
    pointer   operator->() const noexcept { return &(it_->second); }

    const_grid_iterator operator++()    noexcept { it_ = next(it_,level_); return *this; }
    const_grid_iterator operator++(int) noexcept
    {
        const iter old = it_;
        it_ = next(it_,level_);
        return const_grid_iterator { old, level_ };
    }
    
    static iter next( iter it, uchar level ) noexcept
    {
        const multigrid<order>& m { *(it->second.m_) };
        if ( it->second.level_ > 0 )
        {
            // Find sibling.
            const geoid  id     = it->second.id();
            const geoid  parent = it->second.parent_;
            const iter   daddy  = m.cells_[ parent.lvl ].find(parent);
            const uchar  nkids  = daddy->second.num_children();
            const geoid* begin  = daddy->second.children_.data();
            size_t child_number = std::find( begin, begin + nkids, id ) - begin;

            if ( child_number + 1 < nkids )
            {
                // Sibling exists. Traverse down.
                const geoid sibling = daddy->second.children_[ child_number + 1 ];
                it = m.cells_[ sibling.lvl ].find(sibling);
                return leftmost(it,level);
            }
            else
            {
                // Sibling does not exist. Traverse up.
                return next(daddy,level);
            }
        }
        else
        {
            // Top-level. Just increment iterator and traverse down,
            // unless we have reached the end.
            ++it;
            return (it != m.cells_[0].end()) ? leftmost(it,level) : it;
        }
    }

    static iter leftmost( iter it, uchar level ) noexcept
    {
        const multigrid<order>& m { *(it->second.m_) };
        if ( contains(it,level) )
        {
            return it;
        }
        else
        {
            const geoid child = it->second.children_[0];
            return leftmost( m.cells_[ child.lvl ].find(child), level );
        }
    }

    static bool contains( iter it, uchar level ) noexcept
    {
        return ( it->second.level() == level) ||
               ( it->second.level() <= level && it->second.is_leaf() );
    }
};

template <size_t order>
class multigrid<order>::grid_iterator
{
public:
    using iter = typename std::unordered_map<geoid,tetrahedron<order>>::iterator;
    using size_type         = size_t;
    using difference_type   = size_t;
    using value_type        = tetrahedron<order>;
    using reference         = tetrahedron<order>&;
    using pointer           = tetrahedron<order>*;
    using iterator_category = std::forward_iterator_tag;


    iter  it_;
    uchar level_;

    bool operator==( const grid_iterator &rhs ) const noexcept
    { return level_ == rhs.level_ && it_ == rhs.it_; }

    bool operator!=( const grid_iterator &rhs ) const noexcept
    { return level_ != rhs.level_ || it_ != rhs.it_; }

    bool operator==( const const_grid_iterator &rhs ) const noexcept
    { return level_ == rhs.level_ && it_ == rhs.it_; }

    bool operator!=( const const_grid_iterator &rhs ) const noexcept
    { return level_ != rhs.level_ || it_ != rhs.it_; }

    operator const_grid_iterator() const noexcept
    { return const_grid_iterator { it_, level_ }; }

    reference operator* () const noexcept { return   it_->second;  }
    pointer   operator->() const noexcept { return &(it_->second); }

    grid_iterator operator++()    noexcept { it_ = next(it_,level_); return *this; }
    grid_iterator operator++(int) noexcept
    {
        const iter old = it_;
        it_ = next(it_,level_);
        return grid_iterator { old, level_ };
    }
    
    static iter next( iter it, uchar level ) noexcept
    {
        multigrid<order>& m { *(it->second.m_) };
        if ( it->second.level_ > 0 )
        {
            // Find sibling.
            const geoid  id     = it->second.id();
            const geoid  parent = it->second.parent_;
            const iter   daddy  = m.cells_[ parent.lvl ].find(parent);
            const uchar  nkids  = daddy->second.num_children();
            const geoid* begin  = daddy->second.children_.data();
            const size_t child_number = std::find( begin, begin + nkids, id ) - begin;

            if ( child_number + 1 < nkids )
            {
                // Sibling exists. Traverse down.
                const geoid sibling = daddy->second.children_[ child_number + 1 ];
                it = m.cells_[ sibling.lvl ].find(sibling);
                return leftmost(it,level);
            }
            else
            {
                // Sibling does not exist. Traverse up.
                return next(daddy,level);
            }
        }
        else
        {
            // Top-level. Just increment iterator and traverse down,
            // unless we have reached the end.
            ++it;
            return (it != m.cells_[0].end()) ? leftmost(it,level) : it;
        }
    }

    static iter leftmost( iter it, uchar level ) noexcept
    {
        multigrid<order>& m { *(it->second.m_) };
        if ( contains(it,level) )
        {
            return it;
        }
        else
        {
            const geoid child = it->second.children_[0];
            return leftmost( m.cells_[ child.lvl ].find(child), level );
        }
    }

    static bool contains( iter it, uchar level ) noexcept
    {
        return ( it->second.level() == level) ||
               ( it->second.level() <= level && it->second.is_leaf() );
    }
};

template <size_t order> inline
point edge<order>::Chi( real xi ) const noexcept
{
    using shapefcts1d::N;
    using shapefcts1d::apply;
    return apply<order,point>( nodes_, N<order>( xi ) );
}

template <size_t order> inline
point edge<order>::dChi( real xi ) const noexcept
{
    using shapefcts1d::dNdxi;
    using shapefcts1d::apply;
    return apply<order,point>( nodes_, dNdxi<order>( xi ) );
}

template <size_t order> inline
point edge<order>::ChiAff( real xi ) const noexcept
{
    return (1.-xi)*nodes_.front() + xi*nodes_.back();
}

template <size_t order> inline
point edge<order>::get_node( size_t idx ) const noexcept
{
    return nodes_[ idx ];
}

template <size_t order> inline
geoid edge<order>::id() const noexcept
{
    return geoid { 1, level_, (nodes_[0] + nodes_[order])/2 };
}

template <size_t order> inline
point edge<order>::centre() const noexcept
{
    return Chi( 0.5 );
}

template <size_t order> inline
point triangle<order>::Chi( bary2d pos ) const noexcept
{
    using shapefcts2d::N;
    using shapefcts2d::apply;
    return apply<order,point>( nodes_, N<order>( pos ) );
}

template <size_t order> inline
point triangle<order>::dChidxi( bary2d pos ) const noexcept
{
    using shapefcts2d::dNdxi;
    using shapefcts2d::apply;
    return apply<order,point>( nodes_, dNdxi<order>( pos ) );
}

template <size_t order> inline
point triangle<order>::dChideta( bary2d pos ) const noexcept
{
    using shapefcts2d::dNdeta;
    using shapefcts2d::apply;
    return apply<order,point>( nodes_, dNdeta<order>( pos ) );
}

template <size_t order> inline
point triangle<order>::ChiAff( bary2d pos ) const noexcept
{
    return pos.z0 * nodes_[ 0 ] +
           pos.z1 * nodes_[ 1 ] +
           pos.z2 * nodes_[ 2 ];
}

template <size_t order> inline
point triangle<order>::get_node( size_t idx ) const noexcept
{
    return nodes_[ idx ];
}

template <size_t order> inline
point triangle<order>::normal( bary2d pos ) const noexcept
{
    using shapefcts2d::dNdxi;
    using shapefcts2d::dNdeta;
    using shapefcts2d::apply;
    point n = cross_prod( apply<order,point>(nodes_, dNdxi <order>(pos)),
                          apply<order,point>(nodes_, dNdeta<order>(pos)) );
    return n / n.r();
}

template <size_t order> inline
real triangle<order>::surf_elem( bary2d pos ) const noexcept
{
    using shapefcts2d::dNdxi;
    using shapefcts2d::dNdeta;
    using shapefcts2d::apply;
    point n = cross_prod( apply<order,point>(nodes_, dNdxi <order>(pos)),
                          apply<order,point>(nodes_, dNdeta<order>(pos)) );
    return n.r();
}

template <size_t order> inline
point triangle<order>::centre() const noexcept
{
    return Chi( bary2d { 1./3., 1./3., 1./3. } );
}

template <size_t order> inline
bool triangle<order>::is_on_boundary() const noexcept
{
    return ( cells_[0] == NO_GEOID || cells_[1] == NO_GEOID );
}

template <size_t order> inline
geoid triangle<order>::id() const noexcept
{
    return geoid { 2, level_, (nodes_[0]+nodes_[1]+nodes_[2])/3 };
}

template <size_t order> 
void triangle<order>::link( geoid id ) noexcept
{
    if ( id.lvl == level_ )
    {
        if ( cells_[0] == NO_GEOID )
        {
            cells_[0] = id;
        }
        else
        {
            cells_[1] = id;
        }
    }
    else
    {
        if ( cells_[2] == NO_GEOID )
        {
            cells_[2] = id;
        }
        else
        {
            cells_[3] = id;
        }
    }
}

template <size_t order> inline
point tetrahedron<order>::Chi( bary3d p ) const noexcept
{
    using shapefcts3d::N;
    using shapefcts3d::apply;
    return apply<order,point>( nodes_, N<order>(p) );
}

template <size_t order> inline
tensor tetrahedron<order>::dChi( bary3d p ) const noexcept
{
    using shapefcts3d::dNdxi;
    using shapefcts3d::dNdeta;
    using shapefcts3d::dNdzeta;
    using shapefcts3d::apply;

    point dxi   = apply<order,point>( nodes_, dNdxi  <order>(p) );
    point deta  = apply<order,point>( nodes_, dNdeta <order>(p) );
    point dzeta = apply<order,point>( nodes_, dNdzeta<order>(p) );

    return tensor
    {
        dxi.x, deta.x, dzeta.x,
        dxi.y, deta.y, dzeta.y,
        dxi.z, deta.z, dzeta.z
    };
}

template <size_t order> inline
point tetrahedron<order>::ChiAff( bary3d p ) const noexcept
{
    return p.z0 * nodes_[ 0 ] +
           p.z1 * nodes_[ 1 ] +
           p.z2 * nodes_[ 2 ] +
           p.z3 * nodes_[ 3 ];
}

template <size_t order> inline
real tetrahedron<order>::vol_elem( bary3d p ) const noexcept
{
    return std::abs(dChi(p).det());
}

template <size_t order> inline
bary3d tetrahedron<order>::ChiAffInv( point p ) const noexcept
{
    tensor T
    { 
        nodes_[0].x - nodes_[3].x, nodes_[1].x - nodes_[3].x, nodes_[2].x - nodes_[3].x,
        nodes_[0].y - nodes_[3].y, nodes_[1].y - nodes_[3].y, nodes_[2].y - nodes_[3].y,
        nodes_[0].z - nodes_[3].z, nodes_[1].z - nodes_[3].z, nodes_[2].z - nodes_[3].z
    };
    point l = (T.adj() * (p - nodes_[3]))/T.det();
    return bary3d { l.x, l.y, l.z, 1 - l.x - l.y - l.z };
}

template <size_t order> inline
point tetrahedron<order>::get_node( size_t idx ) const noexcept
{
    return nodes_[ idx ];
}

template <size_t order> inline
triangle<order>& tetrahedron<order>::get_face( size_t idx )
{
    return m_->get_face( faces_[idx] );
}

template <size_t order> inline
const triangle<order>& tetrahedron<order>::get_face( size_t idx ) const
{
    return m_->get_face( faces_[idx] );
}


template <size_t order> inline
edge<order>& tetrahedron<order>::get_edge( size_t idx )
{
    return m_->get_edge( edges_[idx] );
}

template <size_t order> inline
const edge<order>& tetrahedron<order>::get_edge( size_t idx ) const 
{
    return m_->get_edge( edges_[idx] );
}

template <size_t order> inline
geoid tetrahedron<order>::get_face_id( size_t idx ) const noexcept
{
    return faces_[idx];
}

template <size_t order> inline
geoid tetrahedron<order>::get_edge_id( size_t idx ) const noexcept
{
    return edges_[idx];
}

template <size_t order> inline
geoid tetrahedron<order>::id() const noexcept
{
    return geoid { 3, level_, (nodes_[0]+nodes_[1]+nodes_[2]+nodes_[3])/4 };
}

template <size_t order> inline
point tetrahedron<order>::centre() const noexcept
{
    return Chi( bary3d { 1./4., 1./4., 1./4., 1./4. } );
}

template <size_t order> inline
point tetrahedron<order>::pos() const noexcept
{
    return centre();
}

template <size_t order> inline
void tetrahedron<order>::bounding_box( point &min, point &max ) const noexcept
{
    constexpr real MAX = std::numeric_limits<real>::max();
    min.x = min.y = min.z =  MAX;
    max.x = max.y = max.z = -MAX;
    for ( point node: nodes_ )
    {
        min.x = std::min( min.x, node.x ); max.x = std::max( max.x, node.x );
        min.y = std::min( min.y, node.y ); max.y = std::max( max.y, node.y );
        min.z = std::min( min.z, node.z ); max.z = std::max( max.z, node.z );
    }
}

template <size_t order> inline
uchar tetrahedron<order>::level() const noexcept
{
    return level_;
}

template <size_t order> inline
bool tetrahedron<order>::is_leaf() const noexcept
{
    return ref_status_ == NoRef;
}

template <size_t order> inline
bool tetrahedron<order>::is_regular() const noexcept
{
    return ( level_ == 0 ) ||
           ( m_->get_cell(parent_).ref_status_ == RegRef );
}

template <size_t order> inline
uchar tetrahedron<order>::num_children() const noexcept
{
    return refinement_children[ ref_status_ ];
}

template <size_t order>
bool tetrahedron<order>::has_child_with_refmark() const noexcept
{
    for ( uchar i = 0; i < num_children(); ++i )
    {
        const tetrahedron<order>& child { m_->get_cell( children_[ i ] ) };
        if ( child.ref_mark_ == Ref )
            return true;
    }
    return false;
}

template <size_t order> inline
bool tetrahedron<order>::has_edge( geoid id ) const noexcept
{
    for ( geoid i: edges_ )
        if ( i == id )
            return true;
    return false;
}

template <size_t order>
bool tetrahedron<order>::has_refined_child_edge() const noexcept
{
    for ( uchar i = 0; i < num_children(); ++i )
    {
        const tetrahedron& child { m_->get_cell( children_[ i ] ) };
        for ( geoid j: child.edges_ )
        {
            if ( ! has_edge(j) && m_->get_edge(j).counter_ > 0 )
                return true;
        }
    }
    return false;
}

template <size_t order>
uchar tetrahedron<order>::get_edge_ref_pattern() const noexcept
{
    constexpr uchar one = 1;
    uchar result { 0 };
    for ( uchar i = 0; i < 6; ++i )
    {
        if ( m_->get_edge(edges_[i]).counter_ > 0 )
            result |= (one << i);
    }
    return result;
}

template <size_t order> inline
void tetrahedron<order>::set_ref() noexcept
{
    if ( is_leaf() ) ref_mark_ = Ref;
}

template <size_t order> inline
void tetrahedron<order>::set_del() noexcept
{
    if ( is_leaf() ) ref_mark_ = Del;
}

template <size_t order> inline
multigrid<order>::multigrid( const multigrid &rhs ):
edges_ { rhs.edges_ }, faces_ { rhs.faces_ }, cells_ { rhs.cells_ }
{
    reset_multigrid_pointers();
}

template <size_t order> inline
multigrid<order>::multigrid( multigrid &&rhs ) noexcept:
edges_ { std::move(rhs.edges_) },
faces_ { std::move(rhs.faces_) },
cells_ { std::move(rhs.cells_) }
{
    reset_multigrid_pointers();
}

template <size_t order> inline
multigrid<order>& multigrid<order>::operator=( const multigrid &rhs )
{
    edges_ = rhs.edges_;
    faces_ = rhs.faces_;
    cells_ = rhs.cells_;
    reset_multigrid_pointers();
    return *this;
}

template <size_t order> inline
multigrid<order>& multigrid<order>::operator=( multigrid &&rhs ) noexcept
{
    edges_ = std::move( rhs.edges_ );
    faces_ = std::move( rhs.faces_ );
    cells_ = std::move( rhs.cells_ );
    reset_multigrid_pointers();
    return *this;
}

template <size_t order>
void multigrid<order>::reset_multigrid_pointers() noexcept
{
    for ( auto &v: edges_ )
    {
        for ( auto &pair: v )
        {
            pair.second.m_ = this;
        }
    }

    for ( auto &v: faces_ )
    {
        for ( auto &pair: v )
        {
            pair.second.m_ = this;
        }
    }

    for ( auto &v: cells_ )
    {
        for ( auto &pair: v )
        {
            pair.second.m_ = this;
        }
    }
}

template <size_t order> inline
bool  multigrid<order>::has_edge( geoid id ) const noexcept
{
    return find_edge( id ) != nullptr;
}

template <size_t order> inline
edge<order>& multigrid<order>::get_edge( geoid id ) noexcept
{
    return edges_[ id.lvl ].find( id )->second;
}

template <size_t order> inline
const edge<order>& multigrid<order>::get_edge( geoid id ) const noexcept
{
    return edges_[ id.lvl ].find( id )->second;
}

template <size_t order> inline
edge<order>* multigrid<order>::find_edge( geoid id ) noexcept
{
    if ( id.lvl < edges_.size() )
    {
        auto iter = edges_[ id.lvl ].find( id );
        if ( iter != edges_[ id.lvl ].end() )
            return &(iter->second);
    }

    return nullptr;
}

template <size_t order> inline
const edge<order>* multigrid<order>::find_edge( geoid id ) const noexcept
{
    if ( id.lvl < edges_.size() )
    {
        auto iter = edges_[ id.lvl ].find( id );
        if ( iter != edges_[ id.lvl ].end() )
            return &(iter->second);
    }

    return nullptr;
}

template <size_t order> inline
bool multigrid<order>::has_face( geoid id ) const noexcept
{
    return find_face(id) != nullptr;
}

template <size_t order> inline
triangle<order>& multigrid<order>::get_face( geoid id ) noexcept
{
    return faces_[ id.lvl ].find( id )->second;
}

template <size_t order> inline
const triangle<order>& multigrid<order>::get_face( geoid id ) const noexcept
{
    return faces_[ id.lvl ].find( id )->second;
}

template <size_t order> inline
triangle<order>* multigrid<order>::find_face( geoid id ) noexcept
{
    if ( id.lvl < faces_.size() )
    {
        auto iter = faces_[ id.lvl ].find( id );
        if ( iter != faces_[ id.lvl ].end() )
            return &(iter->second);
    }

    return nullptr;
}

template <size_t order> inline
const triangle<order>* multigrid<order>::find_face( geoid id ) const noexcept
{
    if ( id.lvl < faces_.size() )
    {
        auto iter = faces_[ id.lvl ].find( id );
        if ( iter != faces_[ id.lvl ].end() )
            return &(iter->second);
    }

    return nullptr;
}

template <size_t order> inline
bool multigrid<order>::has_cell( geoid id ) const noexcept
{
    return find_cell(id) != nullptr;
}

template <size_t order> inline
tetrahedron<order>& multigrid<order>::get_cell( geoid id ) noexcept
{
    return cells_[ id.lvl ].find(id)->second;
}

template <size_t order> inline
const tetrahedron<order>& multigrid<order>::get_cell( geoid id ) const noexcept
{
    return cells_[ id.lvl ].find(id)->second;
}

template <size_t order> inline
tetrahedron<order>* multigrid<order>::find_cell( geoid id ) noexcept
{
    if ( id.lvl < cells_.size() )
    {
        auto iter = cells_[ id.lvl ].find( id );
        if ( iter != cells_[ id.lvl ].end() )
            return &(iter->second);
    }

    return nullptr;
}

template <size_t order> inline
const tetrahedron<order>* multigrid<order>::find_cell( geoid id ) const noexcept
{
    if ( id.lvl < cells_.size() )
    {
        auto iter = cells_[ id.lvl ].find( id );
        if ( iter != cells_[ id.lvl ].end() )
            return &(iter->second);
    }

    return nullptr;
}

template <size_t order> inline
auto multigrid<order>::level_edges_begin( uchar level )       noexcept 
->      level_edge_iterator
{
    return level_edge_iterator { edges_[ level ].begin() };
}

template <size_t order> inline
auto multigrid<order>::level_edges_end  ( uchar level )       noexcept 
->      level_edge_iterator
{
    return level_edge_iterator { edges_[ level ].end() };
}

template <size_t order> inline
auto multigrid<order>::level_edges_begin( uchar level ) const noexcept 
-> const_level_edge_iterator
{
    return const_level_edge_iterator { edges_[ level ].begin() };
}

template <size_t order> inline
auto multigrid<order>::level_edges_end  ( uchar level ) const noexcept 
-> const_level_edge_iterator
{
    return const_level_edge_iterator { edges_[ level ].end() };
}

template <size_t order> inline
auto multigrid<order>::clevel_edges_begin( uchar level ) const noexcept 
-> const_level_edge_iterator
{
    return const_level_edge_iterator { edges_[ level ].begin() };
}

template <size_t order> inline
auto multigrid<order>::clevel_edges_end  ( uchar level ) const noexcept
-> const_level_edge_iterator
{
    return const_level_edge_iterator { edges_[ level ].end() };
}

template <size_t order> inline
auto multigrid<order>::level_faces_begin( uchar level )       noexcept 
->      level_face_iterator
{
    return level_face_iterator { faces_[ level ].begin() };
}

template <size_t order> inline
auto multigrid<order>::level_faces_end  ( uchar level )       noexcept 
->      level_face_iterator
{
    return level_face_iterator { faces_[ level ].end() };
}

template <size_t order> inline
auto multigrid<order>::level_faces_begin( uchar level ) const noexcept 
-> const_level_face_iterator
{
    return const_level_face_iterator { faces_[ level ].begin() };
}

template <size_t order> inline
auto multigrid<order>::level_faces_end  ( uchar level ) const noexcept 
-> const_level_face_iterator
{
    return const_level_face_iterator { faces_[ level ].end() };
}

template <size_t order> inline
auto multigrid<order>::clevel_faces_begin( uchar level ) const noexcept 
-> const_level_face_iterator
{
    return const_level_face_iterator { faces_[ level ].begin() };
}

template <size_t order> inline
auto multigrid<order>::clevel_faces_end  ( uchar level ) const noexcept
-> const_level_face_iterator
{
    return const_level_face_iterator { faces_[ level ].end() };
}

template <size_t order> inline
auto multigrid<order>::level_cells_begin( uchar level )       noexcept 
->      level_cell_iterator
{
    return level_cell_iterator { cells_[ level ].begin() };
}

template <size_t order> inline
auto multigrid<order>::level_cells_end  ( uchar level )       noexcept 
->      level_cell_iterator
{
    return level_cell_iterator { cells_[ level ].end() };
}

template <size_t order> inline
auto multigrid<order>::level_cells_begin( uchar level ) const noexcept 
-> const_level_cell_iterator
{
    return const_level_cell_iterator { cells_[ level ].begin() };
}

template <size_t order> inline
auto multigrid<order>::level_cells_end  ( uchar level ) const noexcept 
-> const_level_cell_iterator
{
    return const_level_cell_iterator { cells_[ level ].end() };
}

template <size_t order> inline
auto multigrid<order>::clevel_cells_begin( uchar level ) const noexcept 
-> const_level_cell_iterator
{
    return const_level_cell_iterator { cells_[ level ].begin() };
}

template <size_t order> inline
auto multigrid<order>::clevel_cells_end  ( uchar level ) const noexcept
-> const_level_cell_iterator
{
    return const_level_cell_iterator { cells_[ level ].end() };
}

template <size_t order> inline
auto multigrid<order>::grid_begin( uchar level ) noexcept
-> grid_iterator
{
    if ( cells_[0].begin() != cells_[0].end() )
    {
        return grid_iterator
        {
            grid_iterator::leftmost(cells_[0].begin(),level),
            level
        };
    }
    else
    {
        return grid_iterator { cells_[0].end(), level };
    }
}

template <size_t order> inline
auto multigrid<order>::grid_end( uchar level ) noexcept
-> grid_iterator
{
    return grid_iterator { cells_[0].end(), level };
}

template <size_t order> inline
auto multigrid<order>::grid_begin( uchar level ) const noexcept
-> const_grid_iterator
{
    if ( cells_[0].begin() != cells_[0].end() )
    {
        return const_grid_iterator
        {
            const_grid_iterator::leftmost(cells_[0].begin(),level),
            level
        };
    }
    else
    {
        return const_grid_iterator { cells_[0].end(), level };
    }
}

template <size_t order> inline
auto multigrid<order>::grid_end( uchar level ) const noexcept
-> const_grid_iterator
{
    return const_grid_iterator { cells_[0].end(), level };
}

template <size_t order> inline
auto multigrid<order>::grid_cbegin( uchar level ) const noexcept
-> const_grid_iterator
{
    if ( cells_[0].begin() != cells_[0].end() )
    {
        return const_grid_iterator
        {
            const_grid_iterator::leftmost(cells_[0].begin(),level),
            level
        };
    }
    else
    {
        return const_grid_iterator { cells_[0].end(), level };
    }
}

template <size_t order> inline
auto multigrid<order>::grid_cend( uchar level ) const noexcept
-> const_grid_iterator
{
    return const_grid_iterator { cells_[0].end(), level };
}

template <size_t order> inline
uchar multigrid<order>::last_level() const noexcept
{
    return cells_.size() - 1;
}

template <size_t order>
void multigrid<order>::adapt()
{
    const size_t J { cells_.size() - 1 };
    for ( size_t k = J+1; k-- > 0; )
    {
        determine_marks( k );
        marks_for_closure( k );
    }  

    for ( size_t k { 0 }; k <= J; ++k )
    {
        if ( cells_[ k ].size() != 0 )
        {
            if ( k > 0 ) marks_for_closure(k);
            if ( k < J ) unrefine(k);
            refine(k);
        }
    } 

    while ( cells_.back().size() == 0 && cells_.size() > 1 )
    {
        cells_.pop_back();
        faces_.pop_back();
        edges_.pop_back();
    }
}

template <size_t order>
void multigrid<order>::marks_for_closure( uchar level )
{
    for ( auto& pair: cells_[ level ] )
    {
        tetrahedron<order>& t = pair.second;

        if ( t.is_regular() && t.ref_mark_ != RegRef )
        {
            uchar R { t.get_edge_ref_pattern() };
            if ( R == NoRef && t.ref_mark_ == Del )
            {
                // Do nothing.
                // NoRef-mark is set in DetermineMarks(level-1)
            }
            else
            {
                t.ref_mark_ = R;
            }
        }
    }
}

template <size_t order>
void multigrid<order>::determine_marks( uchar level )
{
    for ( auto& pair: cells_[ level ] )
    {
        tetrahedron<order>& t = pair.second;

        if ( t.ref_status_ == NoRef )
        {
            if ( t.is_regular() && t.ref_mark_ == Ref )
            {
                t.ref_mark_ = RegRef;
                for ( geoid i: t.edges_ )
                    ++get_edge( i ).counter_;
            }

            if ( level == 0 && t.ref_mark_ == Del )
            {
                t.ref_mark_ = NoRef;
            }
        }
        else if ( t.ref_status_ == RegRef )
        {
            bool all_del = true;
            for ( geoid i: t.children_ )
                all_del &= ( get_cell(i).ref_mark_ == Del );

            if ( all_del )
            {
                t.ref_mark_ = NoRef;
                for ( geoid i: t.edges_ )
                    --get_edge( i ).counter_;
            }

            for ( geoid i: t.children_ )
                if ( get_cell(i).ref_mark_ == Del )
                    get_cell(i).ref_mark_ = NoRef;
        }
        else // t.status_mark_ == IrregRef.
        {
            if ( t.has_child_with_refmark() ||
                 t.has_refined_child_edge() )
            {
                t.ref_mark_ = RegRef;
                
                for ( geoid i: t.edges_ )
                    ++get_edge( i ).counter_;
            } 
            else
            {
                t.ref_mark_ = NoRef;
            }

            for ( uchar i = 0; i < t.num_children(); ++i )
            {
                get_cell(t.children_[i]).ref_mark_ = NoRef;
            }
        }
    }
}

template <size_t order>
void multigrid<order>::unrefine( uchar level )
{
    std::unordered_set<geoid> edge_delete_list;
    std::unordered_set<geoid> face_delete_list;
    std::unordered_set<geoid> cell_delete_list;

    for ( auto& p: cells_[ level + 1 ] )
        cell_delete_list.insert( p.first );

    for ( auto& p: faces_[ level + 1 ] )
        face_delete_list.insert( p.first );

    for ( auto& p: edges_[ level + 1 ] )
        edge_delete_list.insert( p.first );

    for ( auto& pair: cells_[ level ] )
    {
        tetrahedron<order>& t { pair.second };

        if ( t.ref_status_ != NoRef && t.ref_mark_ == t.ref_status_ )
        {
            for ( size_t i { 0 }; i < t.num_children(); ++i )
            {
                tetrahedron<order>& child { get_cell(t.children_[i]) };
                cell_delete_list.erase( child.id() );
                for ( geoid id: child.faces_ )
                    face_delete_list.erase( id );
                for ( geoid id: child.edges_ )
                    edge_delete_list.erase( id );
            }
        }
    }

    for ( geoid id: cell_delete_list ) cells_[ level + 1 ].erase( id );
    for ( geoid id: face_delete_list ) faces_[ level + 1 ].erase( id );
    for ( geoid id: edge_delete_list ) edges_[ level + 1 ].erase( id );
}

template <size_t order>
void multigrid<order>::refine( uchar level )
{
    const size_t J = cells_.size() - 1;
    if ( level == J )
    {
        cells_.push_back( std::unordered_map<geoid,tetrahedron<order>>() );
        faces_.push_back( std::unordered_map<geoid,triangle<order>>() );
        edges_.push_back( std::unordered_map<geoid,edge<order>>() );
    }

    for ( auto& pair: cells_[ level ] )
    {
        tetrahedron<order>& t { pair.second };
        if ( t.ref_mark_ != t.ref_status_ )
        {
            refine_according_to_mark( t );
            t.ref_status_ = t.ref_mark_;
        }
    }
}

template <size_t order>
void multigrid<order>::refine_according_to_mark( tetrahedron<order> &t )
{
    uchar R = t.get_edge_ref_pattern();
    uchar N = refinement_children[R];
    std::array<uchar,8> rule = refinement_rules[R];

    for ( uchar i = 0; i < N; ++i )
    {
        std::array<std::array<uchar,2>,4> child = child_tetrahedra[rule[i]];
        t.children_[i] = make_child( t, child[0], child[1], child[2], child[3] ).id();
    }
}

template <size_t order> inline
edge<order>& multigrid<order>::add( edge<order> e )
{
    return edges_[ e.level_ ][ e.id() ] = std::move(e);
}

template <size_t order> inline
triangle<order>& multigrid<order>::add( triangle<order> f )
{
    return faces_[ f.level_ ][ f.id() ] = std::move(f);
}

template <size_t order> inline
tetrahedron<order>& multigrid<order>::add( tetrahedron<order> t )
{
    return cells_[ t.level_ ][ t.id() ] = std::move(t);
}

template <size_t order> inline
point multigrid<order>::get_pos( const tetrahedron<order> &t, std::array<uchar,2> p ) const noexcept
{
    constexpr uchar edge_map [4][4]
    {
        // (0,0)  (0,1)  (0,2)  (0,3)
        {   255,     0,     1,     3 },

        // (1,0)  (1,1)  (1,2)  (1,3)
        {     0,   255,     2,     4 },

        // (2,0)  (2,1)  (2,2)  (2,3)
        {     1,     2,   255,     5 },

        // (3,0)  (3,1)  (3,2)  (3,3)
        {     3,     4,     5,   255 }
    };

    if ( p[0] == p[1] )
    {
        return t.nodes_[p[0]];
    }
    else
    {
        if ( order == 1 )
        {
            return (t.nodes_[p[0]] + t.nodes_[p[1]])/2;
        }
        else
        {
            return get_edge(t.edges_[edge_map[p[0]][p[1]]]).centre();
        }
    }
}

template <size_t order> inline
bary3d multigrid<order>::get_bary( std::array<uchar,2> p ) const noexcept
{
    using shapefcts3d::lattice;
    constexpr shapefcts3d::coeffs<1,bary3d> lat { lattice<1>() };
    return (lat[p[0]]+lat[p[1]])/2;
}

template <size_t order>
tetrahedron<order>& multigrid<order>::make_child( const tetrahedron<order> &parent,
                                                  std::array<uchar,2> p0,
                                                  std::array<uchar,2> p1,
                                                  std::array<uchar,2> p2,
                                                  std::array<uchar,2> p3 )
{
    using shapefcts3d::lattice;
    using shapefcts3d::coeffs;
    constexpr coeffs<order,bary3d> lat { lattice<order>() };

    point x0  = get_pos(parent,p0); bary3d bx0 = get_bary(p0);
    point x1  = get_pos(parent,p1); bary3d bx1 = get_bary(p1);
    point x2  = get_pos(parent,p2); bary3d bx2 = get_bary(p2);
    point x3  = get_pos(parent,p3); bary3d bx3 = get_bary(p3);

    tetrahedron<order> t;
    t.m_      = this;
    t.level_  = parent.level_ + 1;
    t.parent_ = parent.id();
    t.nodes_[0] = x0;
    t.nodes_[1] = x1;
    t.nodes_[2] = x2;
    t.nodes_[3] = x3;

    triangle<order> &f0 = make_child_face( parent, p1, p2, p3 );
    triangle<order> &f1 = make_child_face( parent, p0, p2, p3 );
    triangle<order> &f2 = make_child_face( parent, p0, p1, p3 );
    triangle<order> &f3 = make_child_face( parent, p0, p1, p2 );

    t.faces_[0] = f0.id(); f0.link(t.id());
    t.faces_[1] = f1.id(); f1.link(t.id());
    t.faces_[2] = f2.id(); f2.link(t.id());
    t.faces_[3] = f3.id(); f3.link(t.id());

    // Get edges from the faces.
    t.edges_[0] = f2.edges_[2]; // (p0,p1).
    t.edges_[1] = f1.edges_[2]; // (p0,p2).
    t.edges_[2] = f0.edges_[2]; // (p1,p2).
    t.edges_[3] = f1.edges_[1]; // (p0,p3).
    t.edges_[4] = f0.edges_[1]; // (p1,p3).
    t.edges_[5] = f0.edges_[0]; // (p2,p3).

    // Set the higher-order DOFs.
    for ( size_t i = 4; i < lat.size(); ++i )
    {
        bary3d b = lat[i];
        t.nodes_[i] = parent.Chi( b.z0*bx0 + b.z1*bx1 + b.z2*bx2 + b.z3*bx3 );
    }

    return add(t);
}

template <size_t order>
triangle<order>& multigrid<order>::make_child_face( const tetrahedron<order> &parent,
                                                    std::array<uchar,2> p0,
                                                    std::array<uchar,2> p1,
                                                    std::array<uchar,2> p2 )
{
    using shapefcts2d::lattice;
    using shapefcts2d::coeffs;
    constexpr coeffs<order,bary2d> lat { lattice<order>() };

    point x0  = get_pos(parent,p0); bary3d bx0 = get_bary(p0);
    point x1  = get_pos(parent,p1); bary3d bx1 = get_bary(p1);
    point x2  = get_pos(parent,p2); bary3d bx2 = get_bary(p2);

    triangle<order> f;
    f.m_     = this;
    f.level_ = parent.level_ + 1;
    f.nodes_[0] = x0;
    f.nodes_[1] = x1;
    f.nodes_[2] = x2;
   
    triangle<order>* fptr = find_face(f.id()); 
    if ( fptr == nullptr )
    {
        f.edges_[0] = make_child_edge(parent,p1,p2);
        f.edges_[1] = make_child_edge(parent,p0,p2);
        f.edges_[2] = make_child_edge(parent,p0,p1);

        // Set the higher-order DOFs.
        for ( size_t k = 3; k < lat.size(); ++k )
        {
            const bary2d b = lat[k];
            f.nodes_[k] = parent.Chi( b.z0*bx0 + b.z1*bx1 + b.z2*bx2 );
        }
        return add(f);
    }
    else return *fptr;
}

template <size_t order>
geoid multigrid<order>::make_child_edge( const tetrahedron<order> &parent,
                                         std::array<uchar,2> p0,
                                         std::array<uchar,2> p1 )
{
    using shapefcts1d::lattice;
    using shapefcts1d::coeffs;
    constexpr coeffs<order,real> lat { lattice<order>() };

    point x0  = get_pos(parent,p0); bary3d bx0 = get_bary(p0);
    point x1  = get_pos(parent,p1); bary3d bx1 = get_bary(p1);
    point pos = (x0+x1)/2;

    // We only return a geoid. We can thus avoid look-ups in the hash-table
    // if the edge is one of the parent's edges.
    if ( pos == parent.edges_[0].pos ) return parent.edges_[0];
    if ( pos == parent.edges_[1].pos ) return parent.edges_[1];
    if ( pos == parent.edges_[2].pos ) return parent.edges_[2];
    if ( pos == parent.edges_[3].pos ) return parent.edges_[3];
    if ( pos == parent.edges_[4].pos ) return parent.edges_[4];
    if ( pos == parent.edges_[5].pos ) return parent.edges_[5];

    edge<order> e;
    e.m_     = this;
    e.level_ = parent.level_ + 1;
    e.nodes_.front() = x0;
    e.nodes_.back () = x1;

    if ( ! has_edge(e.id()) )
    {
        // Set the higher-order DOFs.
        for ( size_t k = 1; k < order; ++k )
        {
            const real x = lat[k];
            e.nodes_[k] = parent.Chi( (1-x)*bx0 + x*bx1 );
        }
        add(e);
    }
    return e.id();
}

}

