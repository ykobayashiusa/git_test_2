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
#ifndef GEOMETRY_MESH_INDEX_H
#define GEOMETRY_MESH_INDEX_H

#include "geometry/point.h"
#include "geometry/box.h"

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <boost/functional/hash.hpp>

namespace geometry
{

template <typename iterator, typename position_getter>
class mesh_index
{
public:
    mesh_index() = delete;
    mesh_index( real h, position_getter pos = position_getter() );

    mesh_index( iterator begin, iterator end, real h,
                position_getter pos = position_getter() );
    
    void   insert( iterator it );
    size_t erase ( iterator it );

    void update();

    real mesh_size() const noexcept;
    void mesh_size( real h );
 
    template <typename output_iterator>
    void query( box b, output_iterator out ) const;

    template <typename predicate, typename output_iterator>
    void query( box b, predicate pred, output_iterator out ) const;

private:
    struct index_t;
    struct  hash_t;

    index_t index( point p )     const noexcept;
    index_t index( iterator it ) const noexcept;

private:
    real            h_;
    position_getter pos_;

    std::unordered_map<index_t,std::vector<iterator>,hash_t> map_;
};

}

#include "geometry/mesh_index.tpp"
#endif

