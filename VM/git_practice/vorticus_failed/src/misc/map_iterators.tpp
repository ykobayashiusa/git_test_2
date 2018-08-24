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

namespace misc
{

template
<
    typename Key,     // map::key_type
    typename T,       // map::mapped_type
    typename Compare, // map::key_compare
    typename Alloc    // map::allocator_type
>
struct const_map_value_iterator
{
    using base = typename std::map<Key,T,Compare,Alloc>::const_iterator;

    using difference_type   = typename std::iterator_traits<base>::difference_type;
    using value_type        = T;
    using       reference   = const T&;
    using       pointer     = const T*;
    using iterator_category = std::bidirectional_iterator_tag;

    base val {};

    bool operator==( const_map_value_iterator rhs ) const noexcept { return val == rhs.val; }
    bool operator!=( const_map_value_iterator rhs ) const noexcept { return val != rhs.val; }

    reference operator* () const noexcept { return   val->second ; }
    pointer   operator->() const noexcept { return &(val->second); }

    const_map_value_iterator  operator++()    noexcept { return const_map_value_iterator { ++val }; }
    const_map_value_iterator  operator--()    noexcept { return const_map_value_iterator { --val }; }
    const_map_value_iterator  operator++(int) noexcept { return const_map_value_iterator { val++ }; }
    const_map_value_iterator  operator--(int) noexcept { return const_map_value_iterator { val-- }; }
};



template
<
    typename Key,     // map::key_type
    typename T,       // map::mapped_type
    typename Compare, // map::key_compare
    typename Alloc    // map::allocator_type
>
struct map_value_iterator
{
    using base = typename std::map<Key,T,Compare,Alloc>::iterator;

    using difference_type   = typename std::iterator_traits<base>::difference_type;
    using value_type        = T;
    using       reference   = T&;
    using       pointer     = T*;
    using iterator_category = std::bidirectional_iterator_tag;

    base val {};

    bool operator==( map_value_iterator rhs ) const noexcept { return val == rhs.val; }
    bool operator!=( map_value_iterator rhs ) const noexcept { return val != rhs.val; }

    reference operator* () const noexcept { return   val->second ; }
    pointer   operator->() const noexcept { return &(val->second); }

    map_value_iterator  operator++()    noexcept { return map_value_iterator { ++val }; }
    map_value_iterator  operator--()    noexcept { return map_value_iterator { --val }; }
    map_value_iterator  operator++(int) noexcept { return map_value_iterator { val++ }; }
    map_value_iterator  operator--(int) noexcept { return map_value_iterator { val-- }; }

    operator const_map_value_iterator<Key,T,Compare,Alloc>() const noexcept
    {   return const_map_value_iterator<Key,T,Compare,Alloc> { val }; }

    bool operator==( const_map_value_iterator<Key,T,Compare,Alloc> rhs ) const noexcept
    { return rhs == (*this); }

    bool operator!=( const_map_value_iterator<Key,T,Compare,Alloc> rhs ) const noexcept
    { return rhs != (*this); }
};

template
<
    typename Key,     // map::key_type
    typename T,       // map::mapped_type
    typename Compare, // map::key_compare
    typename Alloc    // map::allocator_type
> inline
auto make_value_iterator( typename std::unordered_map<Key,T,Compare,Alloc>::iterator it )
-> map_value_iterator<Key,T,Compare,Alloc>
{
    return map_value_iterator<Key,T,Compare,Alloc> { it };
}

template
<
    typename Key,     // map::key_type
    typename T,       // map::mapped_type
    typename Compare, // map::key_compare
    typename Alloc    // map::allocator_type
> inline
auto make_value_iterator( typename std::unordered_map<Key,T,Compare,Alloc>::const_iterator it )
-> const_map_value_iterator<Key,T,Compare,Alloc>
{
    return const_map_value_iterator<Key,T,Compare,Alloc> { it };
}






template
<
    typename Key,    // unordered_map::key_type
    typename T,      // unordered_map::mapped_type
    typename Hash,   // unordered_map::hasher
    typename Pred,   // unordered_map::key_equal
    typename Alloc   // unordered_map::allocator_type
>
struct const_unordered_map_value_iterator
{
    using base = typename std::unordered_map<Key,T,Hash,Pred,Alloc>::const_iterator;

    using difference_type   = typename std::iterator_traits<base>::difference_type;
    using value_type        = T;
    using       reference   = const T&;
    using       pointer     = const T*;
    using iterator_category = std::forward_iterator_tag;

    base val {};

    bool operator==( const_unordered_map_value_iterator rhs ) const noexcept { return val == rhs.val; }
    bool operator!=( const_unordered_map_value_iterator rhs ) const noexcept { return val != rhs.val; }

    reference operator* () const noexcept { return   val->second ; }
    pointer   operator->() const noexcept { return &(val->second); }

    const_unordered_map_value_iterator operator++()    noexcept { return const_unordered_map_value_iterator { ++val }; }
    const_unordered_map_value_iterator operator++(int) noexcept { return const_unordered_map_value_iterator { val++ }; }
};

template
<
    typename Key,    // unordered_map::key_type
    typename T,      // unordered_map::mapped_type
    typename Hash,   // unordered_map::hasher
    typename Pred,   // unordered_map::key_equal
    typename Alloc   // unordered_map::allocator_type
>
struct unordered_map_value_iterator
{
    using base = typename std::unordered_map<Key,T,Hash,Pred,Alloc>::const_iterator;

    using difference_type   = typename std::iterator_traits<base>::difference_type;
    using value_type        = T;
    using       reference   = T&;
    using       pointer     = T*;
    using iterator_category = std::forward_iterator_tag;

    base val {};

    bool operator==( unordered_map_value_iterator rhs ) const noexcept { return val == rhs.val; }
    bool operator!=( unordered_map_value_iterator rhs ) const noexcept { return val != rhs.val; }

    reference operator* () const noexcept { return   val->second ; }
    pointer   operator->() const noexcept { return &(val->second); }

    unordered_map_value_iterator  operator++()    noexcept { return unordered_map_value_iterator { ++val }; }
    unordered_map_value_iterator  operator++(int) noexcept { return unordered_map_value_iterator { val++ }; }

    operator const_unordered_map_value_iterator<Key,T,Hash,Pred,Alloc>() const noexcept
    {   return const_unordered_map_value_iterator<Key,T,Hash,Pred,Alloc> { val }; }

    bool operator==( const_unordered_map_value_iterator<Key,T,Hash,Pred,Alloc> rhs ) const noexcept
    { return rhs == (*this); }

    bool operator!=( const_unordered_map_value_iterator<Key,T,Hash,Pred,Alloc> rhs ) const noexcept
    { return rhs != (*this); }
};

template
<
    typename Key,    // unordered_map::key_type
    typename T,      // unordered_map::mapped_type
    typename Hash,   // unordered_map::hasher
    typename Pred,   // unordered_map::key_equal
    typename Alloc   // unordered_map::allocator_type
> inline
auto make_value_iterator( typename std::unordered_map<Key,T,Hash,Pred,Alloc>::iterator it )
-> unordered_map_value_iterator<Key,T,Hash,Pred,Alloc>
{
    return unordered_map_value_iterator<Key,T,Hash,Pred,Alloc> { it };
}

template
<
    typename Key,    // unordered_map::key_type
    typename T,      // unordered_map::mapped_type
    typename Hash,   // unordered_map::hasher
    typename Pred,   // unordered_map::key_equal
    typename Alloc   // unordered_map::allocator_type
> inline
auto make_value_iterator( typename std::unordered_map<Key,T,Hash,Pred,Alloc>::const_iterator it )
-> const_unordered_map_value_iterator<Key,T,Hash,Pred,Alloc>
{
    return const_unordered_map_value_iterator<Key,T,Hash,Pred,Alloc> { it };
}

}

