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
#ifndef BEM_MULTIGRID_H
#define BEM_MULTIGRID_H

#include <array>
#include <vector>
#include <bitset>
#include <unordered_set>
#include <unordered_map>

#include "bem/geoid.h"
#include "bem/lattice.h"
#include "bem/shape_functions.h"

namespace bem
{

using point = geometry::point;

template <uint order> class edge;
template <uint order> class triangle;
template <uint order> class multigrid;

constexpr std::bitset<4>  NoRef { 0ull };
constexpr std::bitset<4> RegRef { 7ull };
constexpr std::bitset<4>    Ref { 8ull };
constexpr std::bitset<4>    Del { 9ull };

template <uint order>
class edge
{
    friend class triangle<order>;
    friend class multigrid<order>;

public:
    edge() = default;

    point Chi   ( real xi ) const noexcept;
    point ChiAff( real xi ) const noexcept;
    

    bool  has_triangle( size_t idx ) const noexcept;
    const triangle<order>& get_triangle( size_t idx ) const noexcept;
          triangle<order>& get_triangle( size_t idx )       noexcept;

    point get_node( size_t idx ) const noexcept;

    geoid id()     const noexcept;
    point centre() const noexcept;

    edge& get_child( point node ) const noexcept;

private:
    edge( multigrid<order>* m, point p0, point p1 ) noexcept;

private:
    multigrid<order>*   m_ { nullptr };
    bool          refined_ { false };
    unsigned char counter_ { 0 };
    unsigned char   level_ { 0 };

    std::array<point,order+1>    nodes_    {{}};
    std::array<geoid,2>          children_ {{ NO_GEOID, NO_GEOID }};
    std::array<geoid,4>          trias_    {{ NO_GEOID, NO_GEOID, NO_GEOID, NO_GEOID }};
};

template <uint order>
class triangle
{
    friend class edge<order>;
    friend class multigrid<order>;

public:
    triangle() = default;

    point  Chi    ( bary pos ) const noexcept;
    point  ChiAff ( bary pos ) const noexcept;
    point dChidxi ( bary pos ) const noexcept;
    point dChideta( bary pos ) const noexcept;

    point normal   ( bary pos ) const noexcept;
    real  surf_elem( bary pos ) const noexcept;

    point centre() const noexcept;
    point    pos() const noexcept;
    geoid     id() const noexcept;

    void bounding_box( point& min, point& max ) const noexcept;

    real diam() const noexcept;

    unsigned char level() const noexcept;
    bool        is_leaf() const noexcept;
    bool     is_regular() const noexcept;
    size_t num_children() const noexcept;

    bool has_child_with_refmark() const noexcept;
    bool has_refined_child_edge() const noexcept;
    bool has_edge( geoid p_id )   const noexcept;

    std::bitset<4> get_edge_ref_pattern() const noexcept;

    point         get_node( size_t idx ) const noexcept;
    edge<order>&  get_edge( size_t idx ) const noexcept;

    geoid get_edge_id( size_t idx ) const noexcept;

    void set_ref_mark()   noexcept;
    void set_noref_mark() noexcept;
    void set_del_mark()   noexcept;

    std::bitset<4> get_ref_status() const noexcept;
    geoid          get_child_id( size_t idx ) const noexcept;

private:
    triangle( multigrid<order> *m, point p0, point p1, point p2 ) noexcept;
    

private:
    multigrid<order>*  m_ { nullptr };
    unsigned char  level_ { 0 };
    geoid         parent_ = NO_GEOID;
    std::bitset<4>  ref_status_ {};
    std::bitset<4>  ref_mark_   {};

    lattice::nodes<order> nodes_ {{}};
    std::array<geoid,3>   edges_ {{ NO_GEOID, NO_GEOID, NO_GEOID }};

    std::array<geoid,4> children_ {{ NO_GEOID, NO_GEOID, NO_GEOID, NO_GEOID }};
};

template <uint order>
real dist( const triangle<order>& t1, const triangle<order>& t2 );

template <uint order>
class multigrid
{
public:
    class grid;

    multigrid();

          bool             has_triangle( geoid id ) const noexcept;
          triangle<order>& get_triangle( geoid id )       noexcept;
    const triangle<order>& get_triangle( geoid id ) const noexcept;

          bool         has_edge( geoid id ) const noexcept;
          edge<order>& get_edge( geoid id )       noexcept;
    const edge<order>& get_edge( geoid id ) const noexcept;

    void adapt();

    unsigned char last_level() const noexcept;

private:
    void create_edges();
    void link      ( triangle<order>& t, edge<order>& e );
    void link_child( triangle<order>& t, edge<order>& e );
    void split( edge<order>& e );

    void determine_marks  ( unsigned char level );
    void marks_for_closure( unsigned char level );
    void unrefine         ( unsigned char level );
    void   refine         ( unsigned char level );
    void refine_according_to_mark( triangle<order>& t );
    void refine_regularly   ( triangle<order>& t );
    void refine_single_irreg( triangle<order>& t );
    void refine_double_irreg( triangle<order>& t );

private:
    std::vector<std::unordered_map<geoid,edge<order>>>     edges_;
    std::vector<std::unordered_map<geoid,triangle<order>>> trias_;
};



template <uint order>
class multigrid<order>::grid
{
public:
    using size_type = std::vector<geoid>::size_type;

    grid( multigrid& m, unsigned char level );
    
    class       iterator;
    class const_iterator;

    iterator begin() noexcept;
    iterator end()   noexcept;

    const_iterator begin() const noexcept;
    const_iterator end()   const noexcept;

    size_type size() const noexcept;

          multigrid<order>& get_multi_grid()       noexcept { return *m_; }
    const multigrid<order>& get_multi_grid() const noexcept { return *m_; }

    unsigned char get_level() const noexcept { return level_; }

private:
    void dfs( const triangle<order>& t );
    bool contains( const triangle<order>& t ) const noexcept;

private:
    multigrid<order>*      m_ { nullptr };
    unsigned char      level_   { 0 };
    std::vector<geoid> trias_   { 0 };
};

template <uint order>
using grid = typename multigrid<order>::grid;

/////////////////////
// Iterator types. //
/////////////////////

// These are somewhat verbose, but need to be defined before we can define
// the member functions of the previous classes.

template <uint order>
class multigrid<order>::grid::iterator
{
public:

    using size_type       = std::vector<geoid>::size_type;
    using difference_type = std::vector<geoid>::difference_type;

    iterator() = default;

    triangle<order>& operator*() const noexcept
    {
        return m->get_triangle(*iter);
    }

    triangle<order>* operator->() const noexcept
    {
        return &(m->get_triangle(*iter));
    }

    triangle<order>& operator[]( size_type offset ) const noexcept
    {
        return m->get_triangle(*(iter+offset));
    }

    iterator& operator++() noexcept
    {
        ++iter;
        return *this;
    }

    iterator operator++(int) noexcept
    {
        iterator result(*this);
        ++iter;
        return result;
    } 

    iterator& operator--() noexcept
    {
        --iter;
        return *this;
    }

    iterator operator--(int) noexcept
    {
        iterator result(*this);
        --iter;
        return result;
    }

    iterator& operator+=( size_type num ) noexcept
    {
        iter += num;
        return *this;
    }

    iterator& operator-=( size_type num ) noexcept
    {
        iter -= num;
        return *this;
    }

    iterator operator+( size_type num ) const noexcept
    {
        iterator result(*this);
        result += num;
        return result;
    }

    iterator operator-( size_type num ) const noexcept
    {
        iterator result(*this);
        result -= num;
        return result;
    }

    difference_type operator-( const iterator& rhs ) const noexcept
    {
        return iter - rhs.iter;
    }

    bool operator<( const iterator& rhs ) const noexcept
    {
        return iter < rhs.iter;
    }

    bool operator>( const iterator& rhs ) const noexcept
    {
        return iter > rhs.iter;
    }

    bool operator<=( const iterator& rhs ) const noexcept
    {
        return iter <= rhs.iter;
    }

    bool operator>=( const iterator& rhs ) const noexcept
    {
        return iter >= rhs.iter;
    }

    bool operator==( const iterator& rhs ) const noexcept
    {
        return iter == rhs.iter;
    }

    bool operator!=( const iterator& rhs ) const noexcept
    {
        return iter != rhs.iter;
    }

private:
    friend class multigrid<order>::grid;
    friend class multigrid<order>::grid::const_iterator;
    
    iterator( multigrid<order>* mm, std::vector<geoid>::const_iterator iiter ):
    m { mm }, iter { iiter } {}

private:
    multigrid<order>*                  m    { nullptr };
    std::vector<geoid>::const_iterator iter {}; 
};


template <uint order>
class multigrid<order>::grid::const_iterator
{
public:

    using size_type       = std::vector<geoid>::size_type;
    using difference_type = std::vector<geoid>::difference_type;

    const_iterator() = default;

    const_iterator( multigrid::grid::iterator it ):
    m { it.m }, iter { it.iter }
    {}

    const triangle<order>& operator*() const noexcept
    {
        return m->get_triangle(*iter);
    }

    const triangle<order>* operator->() const noexcept
    {
        return &(m->get_triangle(*iter));
    }

    const triangle<order>& operator[]( size_type offset ) const noexcept
    {
        return m->get_triangle(*(iter+offset));
    }

    const_iterator& operator++() noexcept
    {
        ++iter;
        return *this;
    }

    const_iterator operator++(int) noexcept
    {
        const_iterator result(*this);
        ++iter;
        return result;
    }

    const_iterator& operator--() noexcept
    {
        --iter;
        return *this;
    }

    const_iterator operator--(int) noexcept
    {
        const_iterator result(*this);
        --iter;
        return result;
    }

    const_iterator& operator+=( size_type num ) noexcept
    {
        iter += num;
        return *this;
    }

    const_iterator& operator-=( size_type num ) noexcept
    {
        iter -= num;
        return *this;
    }

    const_iterator operator+( size_type num ) const noexcept
    {
        const_iterator result(*this);
        result += num;
        return result;
    }

    const_iterator operator-( size_type num ) const noexcept
    {
        const_iterator result(*this);
        result -= num;
        return result;
    }

    difference_type operator-( const const_iterator& rhs ) const noexcept
    {
        return iter - rhs.iter;
    }

    bool operator<( const const_iterator& rhs ) const noexcept
    {
        return iter < rhs.iter;
    }

    bool operator>( const const_iterator& rhs ) const noexcept
    {
        return iter > rhs.iter;
    }

    bool operator<=( const const_iterator& rhs ) const noexcept
    {
        return iter <= rhs.iter;
    }

    bool operator>=( const const_iterator& rhs ) const noexcept
    {
        return iter >= rhs.iter;
    }

    bool operator==( const const_iterator& rhs ) const noexcept
    {
        return iter == rhs.iter;
    }

    bool operator!=( const const_iterator& rhs ) const noexcept
    {
        return iter != rhs.iter;
    }

private:
    friend class multigrid<order>::grid;
    const_iterator( multigrid<order>* mm, std::vector<geoid>::const_iterator iiter ):
    m { mm }, iter { iiter } {}

private:
    const multigrid*                   m    { nullptr };
    std::vector<geoid>::const_iterator iter {}; 
};

}


#include "bem/multigrid.tpp"
#endif

