/*
 * Copyright (C) 2016 Matthias Kirchhart
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
#ifndef MISC_MAP_ITERATORS_H
#define MISC_MAP_ITERATORS_H

#include <map>
#include <unordered_map>

namespace misc
{

template
<
    typename Key,                                          // map::key_type
    typename T,                                            // map::mapped_type
    typename Compare = std::less<Key>,                     // map::key_compare
    typename Alloc   = std::allocator<std::pair<const Key,T> >  // map::allocator_type
>
struct map_value_iterator;

template
<
    typename Key,                                          // map::key_type
    typename T,                                            // map::mapped_type
    typename Compare = std::less<Key>,                     // map::key_compare
    typename Alloc   = std::allocator<std::pair<const Key,T> >  // map::allocator_type
>
struct const_map_value_iterator;


template
<
    typename Key,                                         // unordered_map::key_type
    typename T,                                           // unordered_map::mapped_type
    typename Hash  = std::hash<Key>,                      // unordered_map::hasher
    typename Pred  = std::equal_to<Key>,                  // unordered_map::key_equal
    typename Alloc = std::allocator< std::pair<const Key,T> >  // unordered_map::allocator_type
>
struct unordered_map_value_iterator;

template
<
    typename Key,                                         // unordered_map::key_type
    typename T,                                           // unordered_map::mapped_type
    typename Hash  = std::hash<Key>,                      // unordered_map::hasher
    typename Pred  = std::equal_to<Key>,                  // unordered_map::key_equal
    typename Alloc = std::allocator< std::pair<const Key,T> >  // unordered_map::allocator_type
>
struct const_unordered_map_value_iterator;




template
<
    typename Key,                                          // map::key_type
    typename T,                                            // map::mapped_type
    typename Compare = std::less<Key>,                     // map::key_compare
    typename Alloc   = std::allocator<std::pair<const Key,T> >  // map::allocator_type
>
auto make_value_iterator( typename std::map<Key,T,Compare,Alloc>::iterator it )
-> map_value_iterator<Key,T,Compare,Alloc>;


template
<
    typename Key,                                          // map::key_type
    typename T,                                            // map::mapped_type
    typename Compare = std::less<Key>,                     // map::key_compare
    typename Alloc   = std::allocator<std::pair<const Key,T> >  // map::allocator_type
>
auto make_value_iterator( typename std::map<Key,T,Compare,Alloc>::const_iterator it )
-> const_map_value_iterator<Key,T,Compare,Alloc>;

template
<
    typename Key,                                         // unordered_map::key_type
    typename T,                                           // unordered_map::mapped_type
    typename Hash  = std::hash<Key>,                      // unordered_map::hasher
    typename Pred  = std::equal_to<Key>,                  // unordered_map::key_equal
    typename Alloc = std::allocator< std::pair<const Key,T> >  // unordered_map::allocator_type
>
auto make_value_iterator( typename std::unordered_map<Key,T,Hash,Pred,Alloc>::iterator it )
-> unordered_map_value_iterator<Key,T,Hash,Pred,Alloc>;

template
<
    typename Key,                                         // unordered_map::key_type
    typename T,                                           // unordered_map::mapped_type
    typename Hash  = std::hash<Key>,                      // unordered_map::hasher
    typename Pred  = std::equal_to<Key>,                  // unordered_map::key_equal
    typename Alloc = std::allocator< std::pair<const Key,T> >  // unordered_map::allocator_type
>
auto make_value_iterator( typename std::unordered_map<Key,T,Hash,Pred,Alloc>::const_iterator it )
-> const_unordered_map_value_iterator<Key,T,Hash,Pred,Alloc>;

}

#include "misc/map_iterators.tpp"
#endif

