/*
 * Copyright (C) 2017 Matthias Kirchhart
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

namespace geometry
{

template <typename iterator, typename position_getter>
mesh_index<iterator,position_getter>::mesh_index( real h, position_getter pos ):
h_ ( h ), pos_ ( pos )
{}

template <typename iterator, typename position_getter>
mesh_index<iterator,position_getter>::mesh_index( iterator begin, iterator end, real h, position_getter pos ):
h_ ( h ), pos_ ( pos )
{
    while ( begin != end ) insert( begin++ );
}

template <typename iterator, typename position_getter>
void mesh_index<iterator,position_getter>::insert( iterator it )
{
    map_[ index(it) ].push_back(it);
}

template <typename iterator, typename position_getter>
size_t mesh_index<iterator,position_getter>::erase( iterator it )
{
    size_t erased { 0 };
    auto box_iter = map_.find( index(it) );
    if ( box_iter != map_.end() )
    {
        std::vector<iterator>& v = box_iter->second;
        auto delpos = std::find( v.begin(), v.end(), it );
        while ( delpos != v.end() )
        {
            delpos = v.erase( delpos );
            erased++;
            delpos = std::find( delpos, v.end(), it );
        }

        if ( v.empty() )
            map_.erase( box_iter );
    }

    return erased;
}

template <typename iterator, typename position_getter>
void mesh_index<iterator,position_getter>::update()
{
    std::vector<iterator> content;
    for ( auto& b: map_ )
    {
        std::vector<iterator>& v = b.second;
        std::copy( v.begin(), v.end(), std::back_inserter( content ) );
    }
    map_.clear();
    for ( iterator it: content ) insert(it);
}

template <typename iterator, typename position_getter>
real mesh_index<iterator,position_getter>::mesh_size() const noexcept
{
    return h_;
}

template <typename iterator, typename position_getter>
void mesh_index<iterator,position_getter>::mesh_size( real h )
{
    h_ = h;
    update();
}

template <typename iterator, typename position_getter>
template <typename output_iterator>
void mesh_index<iterator,position_getter>::query( box b, output_iterator out ) const
{
    auto contains = [b,this]( iterator it ) noexcept -> bool
    {
        const point p = pos_(*it);
        return b.min.x < p.x && p.x < b.max.x &&
               b.min.y < p.y && p.y < b.max.y &&
               b.min.z < p.z && p.z < b.max.z;
    };

    const index_t min_idx { index( b.min ) };
    const index_t max_idx { index( b.max ) };

    for ( long i = min_idx.i; i <= max_idx.i; ++i )
    for ( long j = min_idx.j; j <= max_idx.j; ++j )
    for ( long k = min_idx.k; k <= max_idx.k; ++k )
    {
        auto box_iter = map_.find( index_t { i, j, k } );
        if ( box_iter != map_.end() )
        {
            const std::vector<iterator>& content = box_iter->second;
            std::copy_if( content.begin(), content.end(), out, contains );
        }
    }
}

template <typename iterator, typename position_getter>
template <typename predicate, typename output_iterator>
void mesh_index<iterator,position_getter>::query( box b, predicate pred, output_iterator out ) const
{
    auto contains = [b,&pred,this]( iterator it ) noexcept -> bool
    {
        const point p = pos_(*it);
        return b.min.x < p.x && p.x < b.max.x &&
               b.min.y < p.y && p.y < b.max.y &&
               b.min.z < p.z && p.z < b.max.z &&
               pred(it);
    };

    const index_t min_idx { index( b.min ) };
    const index_t max_idx { index( b.max ) };

    for ( long i = min_idx.i; i <= max_idx.i; ++i )
    for ( long j = min_idx.j; j <= max_idx.j; ++j )
    for ( long k = min_idx.k; k <= max_idx.k; ++k )
    {
        auto box_iter = map_.find( index_t { i, j, k } );
        if ( box_iter != map_.end() )
        {
            const std::vector<iterator>& content = box_iter->second;
            std::copy_if( content.begin(), content.end(), out, contains );
        }
    }
}

template <typename iterator, typename position_getter>
struct mesh_index<iterator,position_getter>::index_t
{
    bool operator==( index_t rhs ) const noexcept
    { return i == rhs.i && j == rhs.j && k == rhs.k; }

    long int i, j, k;
};

template <typename iterator, typename position_getter>
struct mesh_index<iterator,position_getter>::hash_t
{
    size_t operator()( mesh_index::index_t idx ) const noexcept
    {
        using boost::hash_combine;
        std::hash<long int> long_hash;

        size_t result { 0 };
        hash_combine( result, long_hash(idx.i) );
        hash_combine( result, long_hash(idx.j) );
        hash_combine( result, long_hash(idx.k) );

        return result;
    }
};

template <typename iterator, typename position_getter> inline
auto mesh_index<iterator,position_getter>::
index( point p ) const noexcept -> index_t
{
    const long i = std::lround( p.x / h_ );
    const long j = std::lround( p.y / h_ );
    const long k = std::lround( p.z / h_ );
    return index_t { i, j, k };
}

template <typename iterator, typename position_getter> inline
auto mesh_index<iterator,position_getter>::
index( iterator it ) const noexcept -> index_t
{
    return index( pos_(*it) );
}

}

