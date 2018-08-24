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
#ifndef FEM_MULTIGRID_H
#define FEM_MULTIGRID_H

#include <array>
#include <vector>
#include <unordered_set>
#include <unordered_map>

#include "geometry/point.h"
#include "geometry/tensor.h"

#include "fem/geoid.h"
#include "fem/shape_functions.h"
#include "fem/refinement_rules.h"

#include "misc/map_iterators.h"


namespace fem
{

using namespace geometry;

template <size_t order> class edge;
template <size_t order> class triangle;
template <size_t order> class tetrahedron;
template <size_t order> class multigrid;
template <size_t order> class multigrid_builder;

constexpr uchar  NoRef {  0 };
constexpr uchar RegRef { 63 };
constexpr uchar    Ref { 64 };
constexpr uchar    Del { 65 };


template <size_t order>
class edge
{
public:
    static_assert( order >= 1, "Order must not be zero." );

    point  Chi   ( real xi ) const noexcept;
    point dChi   ( real xi ) const noexcept;
    point  ChiAff( real xi ) const noexcept;

    point get_node( size_t idx ) const noexcept;

    geoid id()     const noexcept;
    point centre() const noexcept;

private:
    friend class triangle         <order>;
    friend class tetrahedron      <order>;
    friend class multigrid        <order>;
    friend class multigrid_builder<order>;
    friend class std::pair<const geoid,edge<order>>;

    edge()               = default;
    edge( const edge&  ) = default;
    edge(       edge&& ) = default;
    edge& operator=( const edge&  ) = default;
    edge& operator=(       edge&& ) = default;


private:
    using nodes_t = shapefcts1d::coeffs<order,point>;

    multigrid<order>*   m_        { nullptr };
    uchar               level_    { 0 };
    size_t              counter_  { 0 };
    nodes_t             nodes_    {};
};


template <size_t order>
class triangle
{
public:
    static_assert( order >= 1, "Order must not be zero." );

    point  Chi    ( bary2d p ) const noexcept;
    point  ChiAff ( bary2d p ) const noexcept;
    point dChidxi ( bary2d p ) const noexcept;
    point dChideta( bary2d p ) const noexcept;

    point get_node( size_t idx ) const noexcept;

    point normal   ( bary2d p ) const noexcept;
    real  surf_elem( bary2d p ) const noexcept;

    geoid     id() const noexcept;
    point centre() const noexcept;

    bool is_on_boundary() const noexcept;

private:
    friend class edge             <order>;
    friend class tetrahedron      <order>;
    friend class multigrid        <order>;
    friend class multigrid_builder<order>;
    friend class std::pair<const geoid,triangle<order>>;

    triangle() = default;
    triangle( const triangle&  ) = default;
    triangle(       triangle&& ) = default;
    triangle& operator=( const triangle&   ) = default;
    triangle& operator=(       triangle&&  ) = default;

    void link( geoid id ) noexcept;


private:
    using nodes_t = shapefcts2d::coeffs<order,point>;
    multigrid<order>*    m_       { nullptr };
    uchar                level_   { 0 };
    nodes_t              nodes_   {};
    std::array<geoid,3>  edges_   {};
    std::array<geoid,4>  cells_   {};
};


template <size_t order>
class tetrahedron
{
public:
    static_assert( order >= 1, "Order must not be zero." );

    point   Chi   ( bary3d p ) const noexcept;
    tensor dChi   ( bary3d p ) const noexcept;
    point   ChiAff( bary3d p ) const noexcept;
    real  vol_elem( bary3d p ) const noexcept;

    bary3d ChiAffInv( point p ) const noexcept;

    point                   get_node( size_t idx ) const noexcept;
          triangle<order>&  get_face( size_t idx );
    const triangle<order>&  get_face( size_t idx ) const;
          edge    <order>&  get_edge( size_t idx );
    const edge    <order>&  get_edge( size_t idx ) const;
    geoid             get_face_id( size_t idx ) const noexcept;
    geoid             get_edge_id( size_t idx ) const noexcept;

    geoid     id() const noexcept;
    point centre() const noexcept;
    point    pos() const noexcept;
    void  bounding_box( point &min, point &max ) const noexcept;

    uchar         level() const noexcept;
    bool        is_leaf() const noexcept;
    bool     is_regular() const noexcept;
    uchar  num_children() const noexcept;

    bool has_edge( geoid id ) const noexcept;
    bool has_child_with_refmark() const noexcept;
    bool has_refined_child_edge() const noexcept;
    uchar  get_edge_ref_pattern() const noexcept;

    void set_ref() noexcept;
    void set_del() noexcept;

private:
    friend class edge             <order>;
    friend class triangle         <order>;
    friend class multigrid        <order>;
    friend class multigrid_builder<order>;
    friend class std::pair<const geoid,tetrahedron<order>>;

    tetrahedron() = default;
    tetrahedron( const tetrahedron&  ) = default;
    tetrahedron(       tetrahedron&& ) = default;
    tetrahedron& operator=( const tetrahedron&  ) = default;
    tetrahedron& operator=(       tetrahedron&& ) = default;

private:
    using nodes_t = shapefcts3d::coeffs<order,point>;
    multigrid<order>*   m_          { nullptr };
    uchar               level_      { 0 };
    uchar               ref_mark_   { 0 };
    uchar               ref_status_ { 0 };
    nodes_t             nodes_      {};
    std::array<geoid,4> faces_      {};
    std::array<geoid,6> edges_      {};
    std::array<geoid,8> children_   {};
    geoid               parent_     {};
};


template <size_t order>
class multigrid
{

public:
    static_assert( order >= 1, "Order must not be zero." );
    class grid;

    class       grid_iterator;
    class const_grid_iterator;
    using       level_edge_iterator = misc::      unordered_map_value_iterator<geoid,edge<order>>;
    using       level_face_iterator = misc::      unordered_map_value_iterator<geoid,triangle<order>>;
    using       level_cell_iterator = misc::      unordered_map_value_iterator<geoid,tetrahedron<order>>;
    using const_level_edge_iterator = misc::const_unordered_map_value_iterator<geoid,edge<order>>;
    using const_level_face_iterator = misc::const_unordered_map_value_iterator<geoid,triangle<order>>;
    using const_level_cell_iterator = misc::const_unordered_map_value_iterator<geoid,tetrahedron<order>>;

    multigrid() = default;
    multigrid( const multigrid&  );
    multigrid(       multigrid&& ) noexcept;
    multigrid& operator=( const multigrid& );
    multigrid& operator=(       multigrid&& ) noexcept;

          bool          has_edge( geoid id ) const noexcept;
          edge<order>&  get_edge( geoid id )       noexcept;
    const edge<order>&  get_edge( geoid id ) const noexcept;
          edge<order>* find_edge( geoid id )       noexcept;
    const edge<order>* find_edge( geoid id ) const noexcept;

          bool              has_face( geoid id ) const noexcept;
          triangle<order>&  get_face( geoid id )       noexcept;
    const triangle<order>&  get_face( geoid id ) const noexcept;
          triangle<order>* find_face( geoid id )       noexcept;
    const triangle<order>* find_face( geoid id ) const noexcept;

          bool                 has_cell( geoid id ) const noexcept;
          tetrahedron<order>&  get_cell( geoid id )       noexcept;
    const tetrahedron<order>&  get_cell( geoid id ) const noexcept;
          tetrahedron<order>* find_cell( geoid id )       noexcept;
    const tetrahedron<order>* find_cell( geoid id ) const noexcept;


    /*!
     * \name Iterators through hierarchical surplus.
     *
     * These iterators iterate through all edges, faces, and cells of the
     * hierarchical surplus of a specified level.
     */
    //@{
          level_edge_iterator  level_edges_begin( uchar level )       noexcept;
          level_edge_iterator  level_edges_end  ( uchar level )       noexcept;
    const_level_edge_iterator  level_edges_begin( uchar level ) const noexcept;
    const_level_edge_iterator  level_edges_end  ( uchar level ) const noexcept;
    const_level_edge_iterator clevel_edges_begin( uchar level ) const noexcept;
    const_level_edge_iterator clevel_edges_end  ( uchar level ) const noexcept;

          level_face_iterator  level_faces_begin( uchar level )       noexcept;
          level_face_iterator  level_faces_end  ( uchar level )       noexcept;
    const_level_face_iterator  level_faces_begin( uchar level ) const noexcept;
    const_level_face_iterator  level_faces_end  ( uchar level ) const noexcept;
    const_level_face_iterator clevel_faces_begin( uchar level ) const noexcept;
    const_level_face_iterator clevel_faces_end  ( uchar level ) const noexcept;

          level_cell_iterator  level_cells_begin( uchar level )       noexcept;
          level_cell_iterator  level_cells_end  ( uchar level )       noexcept;
    const_level_cell_iterator  level_cells_begin( uchar level ) const noexcept;
    const_level_cell_iterator  level_cells_end  ( uchar level ) const noexcept;
    const_level_cell_iterator clevel_cells_begin( uchar level ) const noexcept;
    const_level_cell_iterator clevel_cells_end  ( uchar level ) const noexcept;
    //@}


    /*!
     * \name Iterators through entire grids.
     *
     * These iterators iterate through all cells of a grid of a specified level
     * in the hierarchical multigrid. They generally refer to cells on varying
     * hierarchical surplusses.
     */
    //@{
          grid_iterator grid_begin ( uchar level )       noexcept;
          grid_iterator grid_end   ( uchar level )       noexcept;
    const_grid_iterator grid_begin ( uchar level ) const noexcept;
    const_grid_iterator grid_end   ( uchar level ) const noexcept;
    const_grid_iterator grid_cbegin( uchar level ) const noexcept;
    const_grid_iterator grid_cend  ( uchar level ) const noexcept;
    //@}

    uchar last_level() const noexcept;

    void adapt();

private:
    friend class multigrid_builder<order>;

    void determine_marks( uchar level );
    void marks_for_closure( uchar level );
    void unrefine( uchar level );
    void refine( uchar level );
    void refine_according_to_mark( tetrahedron<order>& t );

    edge       <order>& add( edge       <order> e );
    triangle   <order>& add( triangle   <order> f );
    tetrahedron<order>& add( tetrahedron<order> t );

    point  get_pos ( const tetrahedron<order> &parent, std::array<uchar,2> p ) const noexcept;
    bary3d get_bary( std::array<uchar,2> p ) const noexcept;

    tetrahedron<order>&  make_child( const tetrahedron<order> &parent,
                                     std::array<uchar,2> p0, std::array<uchar,2> p1,
                                     std::array<uchar,2> p2, std::array<uchar,2> p3 );
    triangle<order>&     make_child_face( const tetrahedron<order> &parent,
                                          std::array<uchar,2> p0, std::array<uchar,2> p1, std::array<uchar,2> p2 );
    geoid                make_child_edge( const tetrahedron<order> &parent,
                                          std::array<uchar,2> p0, std::array<uchar,2> p1 );

    void reset_multigrid_pointers() noexcept;

private:
    std::vector<std::unordered_map<geoid,edge       <order>>>  edges_ { 1 }; 
    std::vector<std::unordered_map<geoid,triangle   <order>>>  faces_ { 1 };
    std::vector<std::unordered_map<geoid,tetrahedron<order>>>  cells_ { 1 };
};

template <size_t order> using       grid_iterator = typename multigrid<order>::      grid_iterator;
template <size_t order> using const_grid_iterator = typename multigrid<order>::const_grid_iterator;

}

#include "fem/multigrid.tpp"
#endif

