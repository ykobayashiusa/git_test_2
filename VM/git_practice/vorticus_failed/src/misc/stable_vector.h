/*
 * Copyright (C) 2015 Matthias Kirchhart
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
#ifndef MISC_STABLE_VECTOR_H
#define MISC_STABLE_VECTOR_H

#include <vector>
#include <iterator>

namespace misc
{

/*!
 * \brief A vector-like class providing reasonable iterator stability.
 *
 * This class is a replacement for the std::vector class providing iterator
 * stability and contiguous memory.
 *
 * In vortex methods we often need to insert particles to prevent the forming
 * of holes in the particle field. Inserting particles into a standard vector
 * may cause all iterators to the vector to be invalidated. The stable vector
 * class from the Boost library circumvents this problem by giving up contiguos
 * memory. However, to obtain optimal performance and to make better use of
 * processor caches, contiguous memory is of prime importance. This class is
 * thus a compromise between std::vector and boost::stable_vector: it uses
 * std::vector as an underyling container and implements iterators based on
 * indeces, rather than pointers. These iterators are – unlike pointers and
 * references – stable.
 *
 * Currently does not provide: insert, erase, emplace, relational operators.
 */
template <typename T, typename Alloc = std::allocator<T>>
class stable_vector
{
public:

    // Types.
    using      value_type = T;
    using  allocator_type = Alloc;
    using       reference = typename std::vector<T,Alloc>::      reference;
    using const_reference = typename std::vector<T,Alloc>::const_reference;
    using       pointer   = typename std::vector<T,Alloc>::      pointer;
    using const_pointer   = typename std::vector<T,Alloc>::const_pointer;
    using       size_type = typename std::vector<T,Alloc>::size_type;
    using difference_type = typename std::vector<T,Alloc>::difference_type;

    class       iterator;
    class const_iterator;    
    using       reverse_iterator = std::reverse_iterator<      iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

public:

    // Construction and copying.
    explicit stable_vector( const allocator_type& alloc = allocator_type() );
    explicit stable_vector( size_type n );
             stable_vector( size_type n, const value_type& val,
                            const allocator_type& alloc = allocator_type() ); 

    template <typename input_iterator>
    stable_vector( input_iterator begin, input_iterator end,
                   const allocator_type& alloc = allocator_type() );

    stable_vector ( std::initializer_list<value_type> il,
                    const allocator_type& alloc = allocator_type() );

    stable_vector( const stable_vector&  ) = default;
    stable_vector(       stable_vector&& ) = default;
    stable_vector& operator=( const stable_vector&  ) = default;
    stable_vector& operator=(       stable_vector&& ) = default;

    stable_vector( const stable_vector&  rhs,
                   const allocator_type& alloc = allocator_type() );
    stable_vector(       stable_vector&& rhs,
                   const allocator_type& alloc = allocator_type() );

    // Iterators.
          iterator  begin()       noexcept;
    const_iterator  begin() const noexcept;
    const_iterator cbegin() const noexcept;

          iterator    end()       noexcept;
    const_iterator    end() const noexcept;
    const_iterator   cend() const noexcept;

          reverse_iterator  rbegin()       noexcept;
    const_reverse_iterator  rbegin() const noexcept;
    const_reverse_iterator crbegin() const noexcept;

          reverse_iterator  rend()       noexcept;
    const_reverse_iterator  rend() const noexcept;
    const_reverse_iterator crend() const noexcept;

    // Capacity.
    bool         empty() const noexcept;
    size_type     size() const noexcept;
    size_type max_size() const noexcept;
    size_type capacity() const noexcept;

    void reserve( size_type n );
    void shrink_to_fit();

    void resize( size_type n );
    void resize( size_type n, const value_type& v );
    
    // Access.
          reference front();
    const_reference front() const;
          reference back();
    const_reference back() const;
          reference operator[]( size_type n );
    const_reference operator[]( size_type n ) const;
          reference         at( size_type n );
    const_reference         at( size_type n ) const;
          pointer data()       noexcept;
    const_pointer data() const noexcept; 

    // Modifiers.
    template <class input_iterator>
    void assign( input_iterator first, input_iterator last );
    void assign( size_type n, const value_type& val );
    void assign( std::initializer_list<value_type> il );
    void push_back( const value_type&  val );
    void push_back(       value_type&& val );
    void pop_back();
    void swap( stable_vector<T,Alloc> &v );
    void clear() noexcept;
   
    // Allocator.
    allocator_type get_allocator() const noexcept; 
    
private:
    std::vector<T,Alloc> v_;
};

template <typename T, typename Alloc>
void swap( stable_vector<T,Alloc>& x, stable_vector<T,Alloc>& y );

}

#include "misc/stable_vector.tpp"
#endif

