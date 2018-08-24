/*
 * Copyright (C) 2014 Matthias Kirchhart
 *
 * This file is part of vorticus.
 *
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

#pragma GCC push_options
#pragma GCC optimize( "O3" )
namespace bem
{

template <uint order> inline
edge<order>::edge( multigrid<order>* m, point p0, point p1 ) noexcept: m_ { m }
{
    nodes_.front() = p0;
    nodes_.back () = p1;
}

template <uint order> inline
point edge<order>::Chi( real xi ) const noexcept
{
    using shapefcts1d::N;
    using shapefcts1d::apply;
    return apply<order,point>( nodes_, N<order>( xi ) );
}

template <uint order> inline
point edge<order>::ChiAff( real xi ) const noexcept
{
    return (1.-xi)*nodes_.front() + xi*nodes_.back();
}

template <uint order> inline
bool edge<order>::has_triangle( size_t idx ) const noexcept
{
    return m_->has_triangle( trias_[ idx ] );
}

template <uint order> inline
const triangle<order>& edge<order>::get_triangle( size_t idx ) const noexcept
{
    return m_->get_triangle( trias_[ idx ] );
}

template <uint order> inline
triangle<order>& edge<order>::get_triangle( size_t idx ) noexcept
{
    return m_->get_triangle( trias_[ idx ] );
}

template <uint order> inline
point edge<order>::get_node( size_t idx ) const noexcept
{
    return nodes_[ idx ];
}

template <uint order> inline
geoid edge<order>::id() const noexcept
{
    return geoid { 1, level_, centre() };
}

template <uint order> inline
point edge<order>::centre() const noexcept
{
    return Chi( 0.5 );
}

template <uint order> inline
edge<order>& edge<order>::get_child( point node ) const noexcept
{
         if ( node == nodes_.front() ) return m_->get_edge( children_[ 0 ] );
    else if ( node == nodes_.back () ) return m_->get_edge( children_[ 1 ] );
    else return m_->get_edge( NO_GEOID );
}

template <uint order> inline
triangle<order>::triangle( multigrid<order> *m, point p0, point p1, point p2 ) noexcept:
m_ { m }, nodes_ {{ p0, p1, p2 }}
{}

template <uint order> inline
point triangle<order>::Chi( bary pos ) const noexcept
{
    using shapefcts2d::N;
    using shapefcts2d::apply;
    return apply<order,point>( nodes_, N<order>( pos ) );
}

template <uint order> inline
point triangle<order>::dChidxi( bary pos ) const noexcept
{
    using shapefcts2d::dNdxi;
    using shapefcts2d::apply;
    return apply<order,point>( nodes_, dNdxi<order>( pos ) );
}

template <uint order> inline
point triangle<order>::dChideta( bary pos ) const noexcept
{
    using shapefcts2d::dNdeta;
    using shapefcts2d::apply;
    return apply<order,point>( nodes_, dNdeta<order>( pos ) );
}

template <uint order> inline
point triangle<order>::ChiAff( bary pos ) const noexcept
{
    return pos.z0 * nodes_[ 0 ] +
           pos.z1 * nodes_[ 1 ] +
           pos.z2 * nodes_[ 2 ];
}

template <uint order> inline
point triangle<order>::normal( bary pos ) const noexcept
{
    using shapefcts2d::dNdxi;
    using shapefcts2d::dNdeta;
    using shapefcts2d::apply;
    point n = cross_prod( apply<order,point>(nodes_, dNdxi <order>(pos)),
                          apply<order,point>(nodes_, dNdeta<order>(pos)) );
    return n / n.r();
}

template <uint order> inline
real triangle<order>::surf_elem( bary pos ) const noexcept
{
    using shapefcts2d::dNdxi;
    using shapefcts2d::dNdeta;
    using shapefcts2d::apply;
    point n = cross_prod( apply<order,point>(nodes_, dNdxi <order>(pos)),
                          apply<order,point>(nodes_, dNdeta<order>(pos)) );
    return n.r();
}

template <uint order> inline
point triangle<order>::centre() const noexcept
{
    return Chi( bary { 1./3., 1./3., 1./3. } );
}

template <uint order> inline
geoid triangle<order>::id() const noexcept
{
    return geoid { 2, level_, centre() };
}

template <uint order> inline
point triangle<order>::pos() const noexcept
{
    return centre();
}

template <uint order> inline
void triangle<order>::bounding_box( point& min, point& max ) const noexcept
{
    constexpr real MAX = std::numeric_limits<real>::max();
    min = point( MAX, MAX, MAX );
    max = point(-MAX,-MAX,-MAX );
    for ( size_t j = 0; j < nodes_.size(); ++j )
    {
	    min.x = std::min( min.x, nodes_[ j ].x );
	    min.y = std::min( min.y, nodes_[ j ].y );
	    min.z = std::min( min.z, nodes_[ j ].z );

	    max.x = std::max( max.x, nodes_[ j ].x );
	    max.y = std::max( max.y, nodes_[ j ].y );
	    max.z = std::max( max.z, nodes_[ j ].z );
    }
}

template <uint order> inline
real triangle<order>::diam() const noexcept
{
    using std::max;
    return max( (nodes_[1]-nodes_[0]).r(), (nodes_[2]-nodes_[0]).r() );
}

template <uint order>
real dist( const triangle<order>& t1, const triangle<order>& t2 )
{
    using std::min;
    real mindist { std::numeric_limits<real>::max() };
    for ( size_t i { 0 }; i < 3; ++i )
    {
        for ( size_t j { 0 }; j < 3; ++j )
        {
            mindist = min( mindist, (t1.get_node(i) - t2.get_node(j)).r() );
        }
    }

    return mindist;
}

template <uint order> inline
unsigned char triangle<order>::level() const noexcept
{
    return level_;
}

template <uint order> inline
bool triangle<order>::is_leaf() const noexcept
{
    return ref_status_.none();
}

template <uint order> inline
bool triangle<order>::is_regular() const noexcept
{
    return ( level() == 0 ) ||
           ( m_->get_triangle(parent_).ref_status_ == RegRef );
}

template <uint order> inline
size_t triangle<order>::num_children() const noexcept
{
    switch ( ref_status_.count() )
    {
    case 1:   return 2;
    case 2:   return 3;
    case 3:   return 4;
    default:  return 0;
    }
}

template <uint order>
bool triangle<order>::has_child_with_refmark() const noexcept
{
    for ( size_t i { 0 }; i < num_children(); ++i )
    {
        const triangle& child { m_->get_triangle( children_[ i ] ) };
        if ( child.ref_mark_ == Ref )
            return true;
    }
    return false;
}

template <uint order>
bool triangle<order>::has_refined_child_edge() const noexcept
{
    for ( size_t i { 0 }; i < num_children(); ++i )
    {
        const triangle& child { m_->get_triangle( children_[ i ] ) };
        for ( size_t j { 0 }; j < 3; ++j )
        {
            if ( ! has_edge( child.edges_[ j ] ) &&
                   m_->get_edge( child.edges_[ j ] ).counter_ > 0 )
                return true;
        }
    }
    return false;
}

template <uint order>
bool triangle<order>::has_edge( geoid p_id ) const noexcept
{
    if ( p_id == NO_GEOID )
        return false;

    for ( size_t i { 0 }; i < 3; ++i )
    {
        if ( p_id == edges_[ i ] )
            return true;
    }
    return false;
}

template <uint order> inline
std::bitset<4> triangle<order>::get_edge_ref_pattern() const noexcept
{
    std::bitset<4> result;

    if ( m_->get_edge( edges_[ 0 ] ).counter_ ) result.set( 0 ); 
    if ( m_->get_edge( edges_[ 1 ] ).counter_ ) result.set( 1 ); 
    if ( m_->get_edge( edges_[ 2 ] ).counter_ ) result.set( 2 ); 

    return result;
}

template <uint order> inline
point triangle<order>::get_node( size_t idx ) const noexcept
{
    return nodes_[ idx ];
}

template <uint order> inline
edge<order>& triangle<order>::get_edge( size_t idx ) const noexcept
{
    return m_->get_edge( edges_[ idx ] );
}

template <uint order> inline
geoid triangle<order>::get_edge_id( size_t idx ) const noexcept
{
    return edges_[ idx ];
}

template <uint order> inline
void triangle<order>::set_ref_mark() noexcept
{
    if ( is_leaf() ) ref_mark_ = Ref;
}

template <uint order> inline
void triangle<order>::set_noref_mark() noexcept
{
    if ( is_leaf() ) ref_mark_ = NoRef;
}

template <uint order> inline
void triangle<order>::set_del_mark() noexcept
{
    if ( is_leaf() ) ref_mark_ = Del;
}

template <uint order> inline
std::bitset<4> triangle<order>::get_ref_status() const noexcept
{
    return ref_status_;
}

template <uint order> inline
geoid triangle<order>::get_child_id( size_t idx ) const noexcept
{
    return children_[ idx ];
}

template <uint order> inline
bool multigrid<order>::has_triangle( geoid id ) const noexcept
{
    if ( id.lvl >= trias_.size() ) return false;
    return trias_[ id.lvl ].find( id ) != trias_[ id.lvl ].end();
}
template <uint order> inline
triangle<order>& multigrid<order>::get_triangle( geoid id ) noexcept
{
    return trias_[ id.lvl ].find( id )->second;
}

template <uint order> inline
const triangle<order>& multigrid<order>::get_triangle( geoid id ) const noexcept
{
    return trias_[ id.lvl ].find( id )->second;
}

template <uint order> inline
bool multigrid<order>::has_edge( geoid id ) const noexcept
{
    if ( id.lvl >= edges_.size() ) return false;
    return edges_[ id.lvl ].find( id ) != edges_[ id.lvl ].end();
}

template <uint order> inline
edge<order>& multigrid<order>::get_edge( geoid id ) noexcept
{
    return edges_[ id.lvl ].find( id )->second;
}

template <uint order> inline
const edge<order>& multigrid<order>::get_edge( geoid id ) const noexcept
{
    return edges_[ id.lvl ].find( id )->second;
}

template <uint order> inline
unsigned char multigrid<order>::last_level() const noexcept
{
    return trias_.size() - 1;
}

template <uint order>
multigrid<order>::multigrid():
edges_ { 1 }, trias_ { 1 }
{
    constexpr std::array<point,8> nodes
    {
        point {  0.,  0.,  0. },
        point {  1.,  0.,  0. },
        point {  0.,  1.,  0. },
        point {  1.,  1.,  0. },
        point {  0.,  0.,  1. },
        point {  1.,  0.,  1. },
        point {  0.,  1.,  1. },
        point {  1.,  1.,  1. }
    };

    std::array<triangle<order>,12> trias
    {
        triangle<order> { this, nodes[ 0 ], nodes[ 1 ], nodes[ 2 ] }, 
        triangle<order> { this, nodes[ 1 ], nodes[ 3 ], nodes[ 2 ] }, 
        triangle<order> { this, nodes[ 1 ], nodes[ 5 ], nodes[ 3 ] }, 
        triangle<order> { this, nodes[ 5 ], nodes[ 7 ], nodes[ 3 ] }, 
        triangle<order> { this, nodes[ 5 ], nodes[ 4 ], nodes[ 6 ] }, 
        triangle<order> { this, nodes[ 5 ], nodes[ 6 ], nodes[ 7 ] }, 
        triangle<order> { this, nodes[ 0 ], nodes[ 2 ], nodes[ 4 ] }, 
        triangle<order> { this, nodes[ 2 ], nodes[ 6 ], nodes[ 4 ] }, 
        triangle<order> { this, nodes[ 2 ], nodes[ 3 ], nodes[ 6 ] }, 
        triangle<order> { this, nodes[ 3 ], nodes[ 7 ], nodes[ 6 ] }, 
        triangle<order> { this, nodes[ 0 ], nodes[ 4 ], nodes[ 1 ] }, 
        triangle<order> { this, nodes[ 4 ], nodes[ 5 ], nodes[ 1 ] }
    };

    constexpr auto positions = lattice::lattice<order>();

    for ( auto& t: trias )
    {
        for ( size_t k { 3 }; k < t.nodes_.size(); ++k )
        {
            t.nodes_[ k ]  = t.ChiAff( positions[ k ] );
        }
        trias_[ 0 ][ t.id() ] = t;
    }

    create_edges();
}

template <uint order>
void multigrid<order>::create_edges()
{
    edges_.clear();
    edges_.push_back( std::unordered_map<geoid,edge<order>>() );

    std::unordered_map<geoid,geoid> edges_created;

    for ( auto &p: trias_[ 0 ] )
    {
        triangle<order>& t = p.second;

        for ( char i = 0; i < 3; ++i )
        {
            const point p0 { t.nodes_[ i ] };
            const point p1 { t.nodes_[ ( i + 1 ) % 3 ] };

            const geoid id { 1, 0, (p0+p1)/2 };

            if ( edges_created.find( id ) == edges_created.end() )
            {
                edge<order> e { this, p0, p1 };
                for ( size_t k = 1; k < e.nodes_.size() - 1; ++k )
                {
                    real pos { ((real) k)/((real) order) };
                    e.nodes_[ k ]  = e.ChiAff( pos );
                }
                edges_[ 0 ][ e.id() ] = e;
                edges_created[ id ] = e.id();
            }
            
            link( t, edges_[ 0 ][ edges_created[ id ] ] );
        }
    }
}

template <uint order> inline
void multigrid<order>::link( triangle<order>& t, edge<order>& e )
{
    point e0 { e.nodes_.front() };
    point e1 { e.nodes_.back () };

    // Edge 2:
    if      ( t.nodes_[ 0 ] == e0 && t.nodes_[ 1 ] == e1 )
    {
        t.edges_[ 2 ] = e.id();
        e.trias_[ 0 ] = t.id();
    }
    else if ( t.nodes_[ 1 ] == e0 && t.nodes_[ 0 ] == e1 )
    {
        t.edges_[ 2 ] = e.id();
        e.trias_[ 1 ] = t.id();
    }

    // Edge 0:
    else if ( t.nodes_[ 1 ] == e0 && t.nodes_[ 2 ] == e1 )
    {
        t.edges_[ 0 ] = e.id();
        e.trias_[ 0 ] = t.id();
    }
    else if ( t.nodes_[ 2 ] == e0 && t.nodes_[ 1 ] == e1 )
    {
        t.edges_[ 0 ] = e.id();
        e.trias_[ 1 ] = t.id();
    }

    // Edge 1:
    else if ( t.nodes_[ 2 ] == e0 && t.nodes_[ 0 ] == e1 )
    {
        t.edges_[ 1 ] = e.id();
        e.trias_[ 0 ] = t.id();
    }
    else if ( t.nodes_[ 0 ] == e0 && t.nodes_[ 2 ] == e1 )
    {
        t.edges_[ 1 ] = e.id();
        e.trias_[ 1 ] = t.id();
    }
}

template <uint order> inline
void multigrid<order>::link_child( triangle<order>& t, edge<order>& e )
{
    point e0 { e.nodes_.front() };
    point e1 { e.nodes_.back () };
 
    // Edge 2:
    if      ( t.nodes_[ 0 ] == e0 && t.nodes_[ 1 ] == e1 )
    {
        t.edges_[ 2 ] = e.id();
        e.trias_[ 2 ] = t.id();
    }
    else if ( t.nodes_[ 1 ] == e0 && t.nodes_[ 0 ] == e1 )
    {
        t.edges_[ 2 ] = e.id();
        e.trias_[ 3 ] = t.id();
    }

    // Edge 0:
    else if ( t.nodes_[ 1 ] == e0 && t.nodes_[ 2 ] == e1 )
    {
        t.edges_[ 0 ] = e.id();
        e.trias_[ 2 ] = t.id();
    }
    else if ( t.nodes_[ 2 ] == e0 && t.nodes_[ 1 ] == e1 )
    {
        t.edges_[ 0 ] = e.id();
        e.trias_[ 3 ] = t.id();
    }

    // Edge 1:
    else if ( t.nodes_[ 2 ] == e0 && t.nodes_[ 0 ] == e1 )
    {
        t.edges_[ 1 ] = e.id();
        e.trias_[ 2 ] = t.id();
    }
    else if ( t.nodes_[ 0 ] == e0 && t.nodes_[ 2 ] == e1 )
    {
        t.edges_[ 1 ] = e.id();
        e.trias_[ 3 ] = t.id();
    }
}

template <uint order> inline
void multigrid<order>::split( edge<order>& e )
{
    if ( e.refined_ ) return; // Edge already split.

    edge<order> e0 { this, e.nodes_.front(), e.centre() };
    edge<order> e1 { this, e.centre(), e.nodes_.back()  };

    e0.level_  = e1.level_  = e.level_ + 1;

    for ( size_t k { 1 }; k < order; ++k )
    {
        const real tmp { ((real) k)/((real) order) };
        e0.nodes_[ k ] = e.Chi( (     tmp)/2. );
        e1.nodes_[ k ] = e.Chi( (1. + tmp)/2. ); 
    }

    e.refined_ = true;
    e.children_[ 0 ] = e0.id();
    e.children_[ 1 ] = e1.id();

    edges_[ e.level_ + 1 ][ e0.id() ] = e0;
    edges_[ e.level_ + 1 ][ e1.id() ] = e1;
}

template <uint order>
void multigrid<order>::determine_marks( unsigned char level )
{
    for ( auto& pair: trias_[ level ] )
    {
        triangle<order>& t = pair.second;

        if ( t.ref_status_ == NoRef )
        {
            if ( t.is_regular() && t.ref_mark_ == Ref )
            {
                t.ref_mark_ = RegRef;
                ++get_edge( t.edges_[ 0 ] ).counter_;
                ++get_edge( t.edges_[ 1 ] ).counter_;
                ++get_edge( t.edges_[ 2 ] ).counter_;
            }

            if ( level == 0 && t.ref_mark_ == Del )
            {
                t.ref_mark_ = NoRef;
            }
        }
        else if ( t.ref_status_ == RegRef )
        {
            if ( get_triangle( t.children_[ 0 ] ).ref_mark_ == Del &&
                 get_triangle( t.children_[ 1 ] ).ref_mark_ == Del &&
                 get_triangle( t.children_[ 2 ] ).ref_mark_ == Del &&
                 get_triangle( t.children_[ 3 ] ).ref_mark_ == Del )
            {
                t.ref_mark_ = NoRef;
                --get_edge( t.edges_[ 0 ] ).counter_;
                --get_edge( t.edges_[ 1 ] ).counter_;
                --get_edge( t.edges_[ 2 ] ).counter_;
            }

            if ( get_triangle(t.children_[0]).ref_mark_ == Del ) get_triangle(t.children_[0]).ref_mark_ = NoRef;
            if ( get_triangle(t.children_[1]).ref_mark_ == Del ) get_triangle(t.children_[1]).ref_mark_ = NoRef;
            if ( get_triangle(t.children_[2]).ref_mark_ == Del ) get_triangle(t.children_[2]).ref_mark_ = NoRef;
            if ( get_triangle(t.children_[3]).ref_mark_ == Del ) get_triangle(t.children_[3]).ref_mark_ = NoRef;
        }
        else // t.status_mark_ == IrregRef.
        {
            if ( t.has_child_with_refmark() ||
                 t.has_refined_child_edge() )
            {
                t.ref_mark_ = RegRef;
                ++get_edge( t.edges_[ 0 ] ).counter_;
                ++get_edge( t.edges_[ 1 ] ).counter_;
                ++get_edge( t.edges_[ 2 ] ).counter_;
            } 
            else
            {
                t.ref_mark_ = NoRef;
            }

            for ( size_t i { 0 }; i < t.num_children(); ++i )
            {
                get_triangle(t.children_[i]).ref_mark_ = NoRef;
            }
        }
    }
}

template <uint order>
void multigrid<order>::marks_for_closure( unsigned char level )
{
    for ( auto& pair: trias_[ level ] )
    {
        triangle<order>& t = pair.second;

        if ( t.is_regular() && t.ref_mark_ != RegRef )
        {
            std::bitset<4> R { t.get_edge_ref_pattern() };
            if ( R.none() && t.ref_mark_ == Del )
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

template <uint order>
void multigrid<order>::unrefine( unsigned char level )
{
    std::unordered_set<geoid> edge_delete_list;
    std::unordered_set<geoid> tria_delete_list;

    for ( auto& p: trias_[ level + 1 ] )
        tria_delete_list.insert( p.first );

    for ( auto& p: edges_[ level + 1 ] )
        edge_delete_list.insert( p.first );

    for ( auto& pair: trias_[ level ] )
    {
        triangle<order>& t { pair.second };

        if ( t.ref_status_ != NoRef && t.ref_mark_ == t.ref_status_ )
        {
            for ( size_t i { 0 }; i < t.num_children(); ++i )
            {
                triangle<order>& child { get_triangle(t.children_[i]) };
                tria_delete_list.erase( child.id() );
                edge_delete_list.erase( child.edges_[ 0 ] );
                edge_delete_list.erase( child.edges_[ 1 ] );
                edge_delete_list.erase( child.edges_[ 2 ] );
            }
        }
    }

    for ( geoid id: tria_delete_list ) trias_[ level + 1 ].erase( id );
    for ( geoid id: edge_delete_list ) edges_[ level + 1 ].erase( id );
}

template <uint order>
void multigrid<order>::refine( unsigned char level )
{
    const size_t J = trias_.size() - 1;
    if ( level == J )
    {
        trias_.push_back( std::unordered_map<geoid,triangle<order>>() );
        edges_.push_back( std::unordered_map<geoid,edge<order>>() );
    }

    for ( auto& pair: trias_[ level ] )
    {
        triangle<order>& t { pair.second };
        if ( get_edge( t.edges_[ 0 ] ).counter_ ) split( get_edge( t.edges_[ 0 ] ) );
        if ( get_edge( t.edges_[ 1 ] ).counter_ ) split( get_edge( t.edges_[ 1 ] ) );
        if ( get_edge( t.edges_[ 2 ] ).counter_ ) split( get_edge( t.edges_[ 2 ] ) );
    }

    for ( auto& pair: trias_[ level ] )
    {
        triangle<order>& t { pair.second };
        if ( t.ref_mark_ != t.ref_status_ )
        {
            refine_according_to_mark( t );
            t.ref_status_ = t.ref_mark_;
        }
    }
}

template <uint order>
void multigrid<order>::adapt()
{
    const size_t J { trias_.size() - 1 };
    for ( size_t k = J+1; k-- > 0; )
    {
        determine_marks( k );
        marks_for_closure( k );
    }

    for ( size_t k { 0 }; k <= J; ++k )
    {
        if ( trias_[ k ].size() != 0 )
        {
            if ( k > 0 ) marks_for_closure(k);
            if ( k < J ) unrefine(k);
            refine(k);
        }
    } 

    while ( trias_.back().size() == 0 )
    {
        trias_.pop_back();
        edges_.pop_back();
    }
}

template <uint order>
void multigrid<order>::refine_according_to_mark( triangle<order>& t )
{
         if ( t.ref_mark_.count() == 3 ) refine_regularly( t );
    else if ( t.ref_mark_.count() == 2 ) refine_double_irreg( t );
    else if ( t.ref_mark_.count() == 1 ) refine_single_irreg( t );
}

template <uint order>
void multigrid<order>::refine_regularly( triangle<order>& t )
{
    unsigned char level = t.level_ + 1;

    triangle<order> t0 { this, t.nodes_[ 0 ], t.get_edge(2).centre(), t.get_edge(1).centre() };
    triangle<order> t1 { this, t.nodes_[ 1 ], t.get_edge(0).centre(), t.get_edge(2).centre() };
    triangle<order> t2 { this, t.nodes_[ 2 ], t.get_edge(1).centre(), t.get_edge(0).centre() };
    triangle<order> t3 { this, t.get_edge(0).centre(), t.get_edge(1).centre(), t.get_edge(2).centre() };

    constexpr bary t0_0 { 1., 0., 0. }, t0_1 { .5, .5, .0 }, t0_2 { .5, .0, .5 };
    constexpr bary t1_0 { 0., 1., 0. }, t1_1 { .0, .5, .5 }, t1_2 { .5, .5, .0 };
    constexpr bary t2_0 { 0., 0., 1. }, t2_1 { .5, .0, .5 }, t2_2 { .0, .5, .5 };
    constexpr bary t3_0 { 0., .5, .5 }, t3_1 { .5, .0, .5 }, t3_2 { .5, .5, .0 };
    constexpr auto positions = lattice::lattice<order>();

    for ( size_t k { 3 }; k < t.nodes_.size(); ++k )
    {
        const bary pos { positions[ k ] };
        t0.nodes_[ k ] = t.Chi( pos.z0*t0_0 + pos.z1*t0_1 + pos.z2*t0_2 );
        t1.nodes_[ k ] = t.Chi( pos.z0*t1_0 + pos.z1*t1_1 + pos.z2*t1_2 );
        t2.nodes_[ k ] = t.Chi( pos.z0*t2_0 + pos.z1*t2_1 + pos.z2*t2_2 );
        t3.nodes_[ k ] = t.Chi( pos.z0*t3_0 + pos.z1*t3_1 + pos.z2*t3_2 );
    }

    t0.level_  = t1.level_  = t2.level_  = t3.level_  = level;
    t0.parent_ = t1.parent_ = t2.parent_ = t3.parent_ = t.id();

    t.children_[ 0 ] = t0.id();
    t.children_[ 1 ] = t1.id();
    t.children_[ 2 ] = t2.id();
    t.children_[ 3 ] = t3.id();

    trias_[ level ][ t0.id() ] = t0;
    trias_[ level ][ t1.id() ] = t1;
    trias_[ level ][ t2.id() ] = t2;
    trias_[ level ][ t3.id() ] = t3;

    edge<order> e01 { this, get_edge(t.edges_[0]).centre(), get_edge(t.edges_[1]).centre() };
    edge<order> e12 { this, get_edge(t.edges_[1]).centre(), get_edge(t.edges_[2]).centre() };
    edge<order> e20 { this, get_edge(t.edges_[2]).centre(), get_edge(t.edges_[0]).centre() };

    constexpr bary e01_0 { .0, .5, .5 }, e01_1 { .5, .0, .5 };
    constexpr bary e12_0 { .5, .0, .5 }, e12_1 { .5, .5, .0 };
    constexpr bary e20_0 { .5, .5, .0 }, e20_1 { .0, .5, .5 };

    for ( size_t k { 1 }; k < e01.nodes_.size() - 1; ++k )
    {
        real pos { ((real) k)/((real) order) };
        e01.nodes_[ k ] = t.Chi( (1-pos)*e01_0 + pos*e01_1 );
        e12.nodes_[ k ] = t.Chi( (1-pos)*e12_0 + pos*e12_1 );
        e20.nodes_[ k ] = t.Chi( (1-pos)*e20_0 + pos*e20_1 );
    }

    e01.level_ = e12.level_ = e20.level_ = level;
    edges_[ level ][ e01.id() ] = e01;
    edges_[ level ][ e12.id() ] = e12;
    edges_[ level ][ e20.id() ] = e20;

    link( trias_[ level ][ t0.id() ], t.get_edge(1).get_child( t.nodes_[ 0 ] ) ); 
    link( trias_[ level ][ t0.id() ], t.get_edge(2).get_child( t.nodes_[ 0 ] ) ); 
    link( trias_[ level ][ t0.id() ], get_edge( e12.id() ) );
    
    link( trias_[ level ][ t1.id() ], t.get_edge(0).get_child( t.nodes_[ 1 ] ) ); 
    link( trias_[ level ][ t1.id() ], t.get_edge(2).get_child( t.nodes_[ 1 ] ) ); 
    link( trias_[ level ][ t1.id() ], get_edge( e20.id() ) ); 

    link( trias_[ level ][ t2.id() ], t.get_edge(0).get_child( t.nodes_[ 2 ] ) ); 
    link( trias_[ level ][ t2.id() ], t.get_edge(1).get_child( t.nodes_[ 2 ] ) ); 
    link( trias_[ level ][ t2.id() ], get_edge( e01.id() ) ); 

    link( trias_[ level ][ t3.id() ], get_edge( e01.id() ) );
    link( trias_[ level ][ t3.id() ], get_edge( e12.id() ) );
    link( trias_[ level ][ t3.id() ], get_edge( e20.id() ) );
}

template <uint order>
void multigrid<order>::refine_single_irreg( triangle<order>& t )
{
    unsigned char level = t.level_ + 1;

    triangle<order> t0 { this, point {}, point {}, point {} };
    triangle<order> t1 { this, point {}, point {}, point {} };
    edge<order> e { this, point {}, point {} };

    e.level_ = level;
    t0.level_  = t1.level_  = level;
    t0.parent_ = t1.parent_ = t.id();

         if ( t.ref_mark_.test(0) )
    {
        t0.nodes_[ 0 ] = t.nodes_[ 0 ]; t0.nodes_[ 1 ] = t.nodes_[ 1 ]; t0.nodes_[ 2 ] = t.get_edge(0).centre();
        t1.nodes_[ 0 ] = t.nodes_[ 2 ]; t1.nodes_[ 1 ] = t.nodes_[ 0 ]; t1.nodes_[ 2 ] = t.get_edge(0).centre();
        e.nodes_.front() = t.nodes_[ 0 ];
        e.nodes_.back()  = get_edge( t.edges_[0] ).centre();

        constexpr bary t0_0 {  1, .0, .0 }, t0_1 { .0,  1, .0 }, t0_2 { .0, .5, .5 };
        constexpr bary t1_0 { .0, .0,  1 }, t1_1 {  1, .0, .0 }, t1_2 { .0, .5, .5 };
        constexpr bary  e_0 {  1, .0, .0 },  e_1 { .0, .5, .5 };

        constexpr auto positions = lattice::lattice<order>();

        for ( size_t k { 3 }; k < t.nodes_.size(); ++k )
        {
            const bary pos { positions[ k ] };
            t0.nodes_[ k ] = t.Chi( pos.z0*t0_0 + pos.z1*t0_1 + pos.z2*t0_2 );
            t1.nodes_[ k ] = t.Chi( pos.z0*t1_0 + pos.z1*t1_1 + pos.z2*t1_2 );
        }

        for ( size_t k { 1 }; k < e.nodes_.size() - 1; ++k )
        {
            real pos { ((real) k)/((real) order) };
            e.nodes_[ k ] = t.Chi( (1-pos)*e_0 + pos*e_1 );
        }

        link_child( t0, get_edge( t.edges_[ 2 ] ) );
        link( t0, t.get_edge(0).get_child(t.nodes_[1]) );
        link( t0, e );

        link_child( t1, get_edge( t.edges_[ 1 ] ) );
        link( t1, t.get_edge(0).get_child(t.nodes_[2]) );
        link( t1, e );
    }
    else if ( t.ref_mark_.test(1) )
    {
        t0.nodes_[ 0 ] = t.nodes_[ 0 ]; t0.nodes_[ 1 ] = t.nodes_[ 1 ]; t0.nodes_[ 2 ] = t.get_edge(1).centre();
        t1.nodes_[ 0 ] = t.nodes_[ 1 ]; t1.nodes_[ 1 ] = t.nodes_[ 2 ]; t1.nodes_[ 2 ] = t.get_edge(1).centre();
        e.nodes_.front() = t.nodes_[ 1 ];
        e.nodes_.back()  = get_edge( t.edges_[1] ).centre();

        constexpr bary t0_0 {  1, .0, .0 }, t0_1 { .0,  1, .0 }, t0_2 { .5, .0, .5 };
        constexpr bary t1_0 { .0,  1, .0 }, t1_1 { .0, .0,  1 }, t1_2 { .5, .0, .5 };
        constexpr bary  e_0 { .0,  1, .0 },  e_1 { .5, .0, .5 };

        constexpr auto positions = lattice::lattice<order>();

        for ( size_t k { 3 }; k < t.nodes_.size(); ++k )
        {
            const bary pos { positions[ k ] };
            t0.nodes_[ k ] = t.Chi( pos.z0*t0_0 + pos.z1*t0_1 + pos.z2*t0_2 );
            t1.nodes_[ k ] = t.Chi( pos.z0*t1_0 + pos.z1*t1_1 + pos.z2*t1_2 );
        }

        for ( size_t k { 1 }; k < e.nodes_.size() - 1; ++k )
        {
            real pos { ((real) k)/((real) order) };
            e.nodes_[ k ] = t.Chi( (1-pos)*e_0 + pos*e_1 );
        }

        link_child( t0, t.get_edge(2) );
        link( t0, t.get_edge(1).get_child(t.nodes_[0]) );
        link( t0, e );

        link_child( t1, t.get_edge(0) );
        link( t1, t.get_edge(1).get_child(t.nodes_[2]) );
        link( t1, e );
    }
    else // ( t.ref_mark_.test(2) )
    {
        t0.nodes_[ 0 ] = t.nodes_[ 2 ]; t0.nodes_[ 1 ] = t.nodes_[ 0 ]; t0.nodes_[ 2 ] = t.get_edge(2).centre();
        t1.nodes_[ 0 ] = t.nodes_[ 1 ]; t1.nodes_[ 1 ] = t.nodes_[ 2 ]; t1.nodes_[ 2 ] = t.get_edge(2).centre();
        e.nodes_.front() = t.nodes_[ 2 ];
        e.nodes_.back()  = get_edge( t.edges_[2] ).centre();

        constexpr bary t0_0 { .0, .0,  1 }, t0_1 {  1, .0, .0 }, t0_2 { .5, .5, .0 };
        constexpr bary t1_0 { .0,  1, .0 }, t1_1 { .0, .0,  1 }, t1_2 { .5, .5, .0 };
        constexpr bary  e_0 { .0, 0.,  1 },  e_1 { .5, .5, .0 };

        constexpr auto positions = lattice::lattice<order>();

        for ( size_t k { 3 }; k < t.nodes_.size(); ++k )
        {
            const bary pos { positions[ k ] };
            t0.nodes_[ k ] = t.Chi( pos.z0*t0_0 + pos.z1*t0_1 + pos.z2*t0_2 );
            t1.nodes_[ k ] = t.Chi( pos.z0*t1_0 + pos.z1*t1_1 + pos.z2*t1_2 );
        }

        for ( size_t k { 1 }; k < e.nodes_.size() - 1; ++k )
        {
            real pos { ((real) k)/((real) order) };
            e.nodes_[ k ] = t.Chi( (1-pos)*e_0 + pos*e_1 );
        }

        link_child( t0, t.get_edge(1) );
        link( t0, t.get_edge(2).get_child(t.nodes_[0]) );
        link( t0, e );

        link_child( t1, t.get_edge(0) );
        link( t1, t.get_edge(2).get_child(t.nodes_[1]) );
        link( t1, e );
    }

    trias_[ level ][ t0.id() ] = t0;
    trias_[ level ][ t1.id() ] = t1;
    edges_[ level ][  e.id() ] = e;

    t.children_[ 0 ] = t0.id();
    t.children_[ 1 ] = t1.id();
    t.children_[ 2 ] = NO_GEOID;
    t.children_[ 3 ] = NO_GEOID;

}

template <uint order>
void multigrid<order>::refine_double_irreg( triangle<order>& t )
{
    unsigned char level = t.level_ + 1;

    triangle<order> t0 { this, point {}, point {}, point {} };
    triangle<order> t1 { this, point {}, point {}, point {} };
    triangle<order> tc { this, point {}, point {}, point {} };
    edge<order> e0 { this, point {}, point {} };
    edge<order> e1 { this, point {}, point {} };

    e1.level_ = e0.level_ = level;
    t0.level_  = t1.level_  = tc.level_  = level;
    t0.parent_ = t1.parent_ = tc.parent_ = t.id();

         if ( ! t.ref_mark_.test(2) ) // 011
    {
        t0.nodes_[ 0 ] = t.nodes_[ 0 ]; t0.nodes_[ 1 ] = t.nodes_[ 1 ];          t0.nodes_[ 2 ] = t.get_edge(0).centre();
        tc.nodes_[ 0 ] = t.nodes_[ 0 ]; tc.nodes_[ 1 ] = t.get_edge(0).centre(); tc.nodes_[ 2 ] = t.get_edge(1).centre();
        t1.nodes_[ 0 ] = t.nodes_[ 2 ]; t1.nodes_[ 1 ] = t.get_edge(1).centre(); t1.nodes_[ 2 ] = t.get_edge(0).centre();

        e0.nodes_.front() = t.nodes_[ 0 ];          e0.nodes_.back() = t.get_edge(0).centre();
        e1.nodes_.front() = t.get_edge(0).centre(); e1.nodes_.back() = t.get_edge(1).centre();

        constexpr bary t0_0 {  1, .0, .0 }, t0_1 { .0,  1, .0 }, t0_2 { .0, .5, .5 };
        constexpr bary tc_0 {  1, .0, .0 }, tc_1 { .0, .5, .5 }, tc_2 { .5, .0, .5 };
        constexpr bary t1_0 { .0, .0,  1 }, t1_1 { .5, .0, .5 }, t1_2 { .0, .5, .5 };
        constexpr bary e0_0 {  1, .0, .0 }, e0_1 { .0, .5, .5 };
        constexpr bary e1_0 { .0, .5, .5 }, e1_1 { .5, .0, .5 };

        constexpr auto positions = lattice::lattice<order>();

        for ( size_t k { 3 }; k < t.nodes_.size(); ++k )
        {
            const bary pos { positions[ k ] };
            t0.nodes_[ k ] = t.Chi( pos.z0*t0_0 + pos.z1*t0_1 + pos.z2*t0_2 );
            tc.nodes_[ k ] = t.Chi( pos.z0*tc_0 + pos.z1*tc_1 + pos.z2*tc_2 );
            t1.nodes_[ k ] = t.Chi( pos.z0*t1_0 + pos.z1*t1_1 + pos.z2*t1_2 );
        }

        for ( size_t k { 1 }; k < e0.nodes_.size() - 1; ++k )
        {
            real pos { ((real) k)/((real) order) };
            e0.nodes_[ k ] = t.Chi( (1-pos)*e0_0 + pos*e0_1 );
            e1.nodes_[ k ] = t.Chi( (1-pos)*e1_0 + pos*e1_1 );
        }

        link_child( t0, t.get_edge(2) );
        link( t0, t.get_edge(0).get_child(t.nodes_[1]) );
        link( t0, e0 );

        link( tc, e0 );
        link( tc, e1 );
        link( tc, t.get_edge(1).get_child(t.nodes_[0]) );

        link( t1, t.get_edge(0).get_child(t.nodes_[2] ) );
        link( t1, t.get_edge(1).get_child(t.nodes_[2] ) );
        link( t1, e1 );
    }
    else if ( ! t.ref_mark_.test(1) ) // 101
    {
        t0.nodes_[ 0 ] = t.nodes_[ 1 ]; t0.nodes_[ 1 ] = t.get_edge(0).centre(); t0.nodes_[ 2 ] = t.get_edge(2).centre();
        tc.nodes_[ 0 ] = t.nodes_[ 0 ]; tc.nodes_[ 1 ] = t.get_edge(2).centre(); tc.nodes_[ 2 ] = t.get_edge(0).centre();
        t1.nodes_[ 0 ] = t.nodes_[ 2 ]; t1.nodes_[ 1 ] = t.nodes_[ 0 ];          t1.nodes_[ 2 ] = t.get_edge(0).centre();
        e0.nodes_.front() = t.get_edge(2).centre(); e0.nodes_.back() = t.get_edge(0).centre();
        e1.nodes_.front() = t.nodes_[ 0 ];          e1.nodes_.back() = t.get_edge(0).centre();

        constexpr bary t0_0 { .0,  1, .0 }, t0_1 { .0, .5, .5 }, t0_2 { .5, .5, .0 };
        constexpr bary tc_0 {  1, .0, .0 }, tc_1 { .5, .5, .0 }, tc_2 { .0, .5, .5 };
        constexpr bary t1_0 { .0, .0,  1 }, t1_1 {  1, .0, .0 }, t1_2 { .0, .5, .5 };
        constexpr bary e0_0 { .5, .5, .0 }, e0_1 { .0, .5, .5 };
        constexpr bary e1_0 {  1, .0, .0 }, e1_1 { .0, .5, .5 };

        constexpr auto positions = lattice::lattice<order>();

        for ( size_t k { 3 }; k < t.nodes_.size(); ++k )
        {
            const bary pos { positions[ k ] };
            t0.nodes_[ k ] = t.Chi( pos.z0*t0_0 + pos.z1*t0_1 + pos.z2*t0_2 );
            tc.nodes_[ k ] = t.Chi( pos.z0*tc_0 + pos.z1*tc_1 + pos.z2*tc_2 );
            t1.nodes_[ k ] = t.Chi( pos.z0*t1_0 + pos.z1*t1_1 + pos.z2*t1_2 );
        }

        for ( size_t k { 1 }; k < e0.nodes_.size() - 1; ++k )
        {
            real pos { ((real) k)/((real) order) };
            e0.nodes_[ k ] = t.Chi( (1-pos)*e0_0 + pos*e0_1 );
            e1.nodes_[ k ] = t.Chi( (1-pos)*e1_0 + pos*e1_1 );
        }

        link( t0, e0 );
        link( t0, t.get_edge(2).get_child(t.nodes_[1]) );
        link( t0, t.get_edge(0).get_child(t.nodes_[1]) );

        link( tc, e0 );
        link( tc, e1 );
        link( tc, t.get_edge(2).get_child(t.nodes_[0]) );

        link_child( t1, t.get_edge(1) );
        link( t1, t.get_edge(0).get_child(t.nodes_[2]) );
        link( t1, e1 );
    }
    else // 110
    {
        t0.nodes_[ 0 ] = t.nodes_[ 0 ]; t0.nodes_[ 1 ] = t.get_edge(2).centre(); t0.nodes_[ 2 ] = t.get_edge(1).centre();
        tc.nodes_[ 0 ] = t.nodes_[ 1 ]; tc.nodes_[ 1 ] = t.get_edge(1).centre(); tc.nodes_[ 2 ] = t.get_edge(2).centre();
        t1.nodes_[ 0 ] = t.nodes_[ 2 ]; t1.nodes_[ 1 ] = t.get_edge(1).centre(); t1.nodes_[ 2 ] = t.nodes_[ 1 ];
        e0.nodes_.front() = t.get_edge(1).centre(); e0.nodes_.back() = t.get_edge(2).centre();
        e1.nodes_.front() = t.get_edge(1).centre(); e1.nodes_.back() = t.nodes_[ 1 ];

        constexpr bary t0_0 {  1, .0, .0 }, t0_1 { .5, .5, .0 }, t0_2 { .5, .0, .5 };
        constexpr bary tc_0 { .0,  1, .0 }, tc_1 { .5, .0, .5 }, tc_2 { .5, .5, .0 };
        constexpr bary t1_0 { .0, .0,  1 }, t1_1 { .5, .0, .5 }, t1_2 { .0,  1, .0 };
        constexpr bary e0_0 { .5, .0, .5 }, e0_1 { .5, .5, .0 };
        constexpr bary e1_0 { .5, .0, .5 }, e1_1 { .0,  1, .0 };

        constexpr auto positions = lattice::lattice<order>();

        for ( size_t k { 3 }; k < t.nodes_.size(); ++k )
        {
            const bary pos { positions[ k ] };
            t0.nodes_[ k ] = t.Chi( pos.z0*t0_0 + pos.z1*t0_1 + pos.z2*t0_2 );
            tc.nodes_[ k ] = t.Chi( pos.z0*tc_0 + pos.z1*tc_1 + pos.z2*tc_2 );
            t1.nodes_[ k ] = t.Chi( pos.z0*t1_0 + pos.z1*t1_1 + pos.z2*t1_2 );
        }

        for ( size_t k { 1 }; k < e0.nodes_.size() - 1; ++k )
        {
            real pos { ((real) k)/((real) order) };
            e0.nodes_[ k ] = t.Chi( (1-pos)*e0_0 + pos*e0_1 );
            e1.nodes_[ k ] = t.Chi( (1-pos)*e1_0 + pos*e1_1 );
        }

        link( t0, e0 );
        link( t0, t.get_edge(2).get_child( t.nodes_[0] ) );
        link( t0, t.get_edge(1).get_child( t.nodes_[0] ) );

        link( tc, e0 );
        link( tc, e1 );
        link( tc, t.get_edge(2).get_child( t.nodes_[1] ) );

        link( t1, e1 );
        link( t1, t.get_edge(1).get_child( t.nodes_[2] ) );
        link_child( t1, t.get_edge(0) );
    }

    t.children_[ 0 ] = t0.id();
    t.children_[ 1 ] = tc.id();
    t.children_[ 2 ] = t1.id();
    t.children_[ 3 ] = NO_GEOID;

    trias_[ level ][ t0.id() ] = t0;
    trias_[ level ][ tc.id() ] = tc;
    trias_[ level ][ t1.id() ] = t1;
    edges_[ level ][ e0.id() ] = e0;
    edges_[ level ][ e1.id() ] = e1;
}


template <uint order>
multigrid<order>::grid::grid( multigrid<order> &m, unsigned char level ):
m_ { &m }, level_ { level }, trias_ { 0 }
{
    if ( m_->trias_.size() > 0 )
    {
        for ( auto& pair: m_->trias_[ 0 ] )
        {
            dfs( pair.second );
        }
    }
}

template <uint order>
void multigrid<order>::grid::dfs( const triangle<order> &t )
{
    if ( contains(t) )
    {
        trias_.push_back( t.id() );
        return;
    }
    else
    {
        for ( size_t i { 0 }; i < t.num_children(); ++i )
        {
            dfs( m_->get_triangle( t.children_[i] ) );
        }
    }
}

template <uint order> inline
bool multigrid<order>::grid::contains( const triangle<order>& t ) const noexcept
{
    return    ( t.level() == level_ ) 
           || ( t.level() <  level_ && t.is_leaf() );
}

template <uint order> inline
auto multigrid<order>::grid::begin() noexcept -> iterator
{
    return iterator( m_, trias_.begin() );
}

template <uint order> inline
auto multigrid<order>::grid::begin() const noexcept -> const_iterator
{
    return const_iterator( m_, trias_.begin() );
}

template <uint order> inline
auto multigrid<order>::grid::end() noexcept -> iterator
{
    return iterator( m_, trias_.end() );
}

template <uint order> inline
auto multigrid<order>::grid::end() const noexcept -> const_iterator
{
    return const_iterator( m_, trias_.end() );
}

template <uint order> inline
auto multigrid<order>::grid::size() const noexcept -> size_type
{
    return trias_.size();
}

}
#pragma GCC pop_options

