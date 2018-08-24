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

namespace misc
{

template <typename T, typename Alloc>
class stable_vector<T,Alloc>::const_iterator
{
public:
    using value_type        = T;
    using size_type         = typename stable_vector<T,Alloc>::size_type;
    using difference_type   = typename stable_vector<T,Alloc>::difference_type;
    using reference         = typename stable_vector<T,Alloc>::const_reference;
    using pointer           = typename stable_vector<T,Alloc>::const_pointer;
    using iterator_category = std::random_access_iterator_tag;

    friend class stable_vector<T,Alloc>;
    friend class stable_vector<T,Alloc>::iterator;

    const_iterator() noexcept = default;
    const_iterator( const const_iterator&  rhs ) noexcept = default;
    const_iterator(       const_iterator&& rhs ) noexcept = default;
    ~const_iterator() = default;

    const_iterator& operator=( const const_iterator&  rhs ) noexcept = default;
    const_iterator& operator=(       const_iterator&& rhs ) noexcept = default;

    bool operator==( const_iterator rhs ) const noexcept
    { return v == rhs.v && pos == rhs.pos; }

    bool operator!=( const_iterator rhs ) const noexcept
    { return !((*this) == rhs); }

    const_iterator& operator++()    noexcept { ++pos; return *this; }
    const_iterator& operator--()    noexcept { --pos; return *this; }
    const_iterator  operator++(int) noexcept { const_iterator old(*this); ++pos; return old; }
    const_iterator  operator--(int) noexcept { const_iterator old(*this); --pos; return old; }

    const_iterator& operator+=( difference_type n ) noexcept { pos += n; return *this; }
    const_iterator& operator-=( difference_type n ) noexcept { pos -= n; return *this; }

    const_iterator operator+( difference_type n ) const noexcept { const_iterator res(*this); res += n; return res; }
    const_iterator operator-( difference_type n ) const noexcept { const_iterator res(*this); res -= n; return res; }

    difference_type operator-( const_iterator rhs ) const noexcept
    { difference_type res = pos; res -= rhs.pos; return res; }

    reference operator[]( difference_type n ) const noexcept { return (*v)[ pos + n ]; }

    bool operator< ( const_iterator rhs ) const noexcept { return pos <  rhs.pos; }
    bool operator> ( const_iterator rhs ) const noexcept { return pos >  rhs.pos; }
    bool operator<=( const_iterator rhs ) const noexcept { return pos <= rhs.pos; }
    bool operator>=( const_iterator rhs ) const noexcept { return pos >= rhs.pos; }

    reference operator* () const noexcept { return   (*v)[pos]; }
    pointer   operator->() const noexcept { return &((*v)[pos]); }

private:
    const_iterator( const std::vector<T,Alloc> *vv, size_type p ) noexcept: v { vv }, pos { p } {}

    const std::vector<T,Alloc>* v   { nullptr };
    size_type                   pos { 0 };
};

template <typename T, typename Alloc> inline
typename stable_vector<T,Alloc>::const_iterator 
operator+( typename stable_vector<T,Alloc>::difference_type n,
           typename stable_vector<T,Alloc>::const_iterator it ) noexcept
{
    it += n;
    return it;
}




template <typename T, typename Alloc>
class stable_vector<T,Alloc>::iterator
{
public:
    using value_type        = T;
    using size_type         = typename stable_vector<T,Alloc>::size_type;
    using difference_type   = typename stable_vector<T,Alloc>::difference_type;
    using reference         = typename stable_vector<T,Alloc>::reference;
    using pointer           = typename stable_vector<T,Alloc>::pointer;
    using iterator_category = std::random_access_iterator_tag;

    using const_iterator = typename stable_vector<T,Alloc>::const_iterator;

    friend class stable_vector<T,Alloc>;

    iterator() noexcept = default;
    iterator( const iterator&  rhs ) noexcept = default;
    iterator(       iterator&& rhs ) noexcept = default;
    ~iterator() = default;

    iterator& operator=( const iterator&  rhs ) noexcept = default;
    iterator& operator=(       iterator&& rhs ) noexcept = default;

    bool operator==( iterator rhs ) const noexcept
    { return v == rhs.v && pos == rhs.pos; }

    bool operator!=( iterator rhs ) const noexcept
    { return !((*this) == rhs); }

    iterator& operator++()    noexcept { ++pos; return *this; }
    iterator& operator--()    noexcept { --pos; return *this; }
    iterator  operator++(int) noexcept { iterator old(*this); ++pos; return old; }
    iterator  operator--(int) noexcept { iterator old(*this); --pos; return old; }

    iterator& operator+=( difference_type n ) noexcept { pos += n; return *this; }
    iterator& operator-=( difference_type n ) noexcept { pos -= n; return *this; }

    iterator operator+( difference_type n ) const noexcept { iterator res(*this); res += n; return res; }
    iterator operator-( difference_type n ) const noexcept { iterator res(*this); res -= n; return res; }

    difference_type operator-( iterator rhs ) const noexcept
    { difference_type res = pos; res -= rhs.pos; return res; }

    reference operator[]( difference_type n ) const noexcept { return (*v)[ pos + n ]; }

    bool operator< ( iterator rhs ) const noexcept { return pos <  rhs.pos; }
    bool operator> ( iterator rhs ) const noexcept { return pos >  rhs.pos; }
    bool operator<=( iterator rhs ) const noexcept { return pos <= rhs.pos; }
    bool operator>=( iterator rhs ) const noexcept { return pos >= rhs.pos; }

    reference operator* () const noexcept { return   (*v)[pos]; }
    pointer   operator->() const noexcept { return &((*v)[pos]); }

    operator const_iterator() const noexcept { return const_iterator(v,pos); }

private:
    iterator( std::vector<T,Alloc> *vv, size_type p ) noexcept: v { vv }, pos { p } {}

    std::vector<T,Alloc>* v   { nullptr };
    size_type             pos { 0 };
};

template <typename T, typename Alloc> inline
typename stable_vector<T,Alloc>::iterator 
operator+( typename stable_vector<T,Alloc>::difference_type n,
           typename stable_vector<T,Alloc>::iterator it ) noexcept
{
    it += n;
    return it;
}




template <typename T, typename Alloc>
stable_vector<T,Alloc>::stable_vector( const allocator_type& alloc ):
v_(alloc)
{}

template <typename T, typename Alloc>
stable_vector<T,Alloc>::stable_vector( size_type n ):
v_(n)
{}

template <typename T, typename Alloc>
stable_vector<T,Alloc>::stable_vector( size_type n, const value_type& val,
                                       const allocator_type& alloc ):
v_(n,val,alloc)
{}

template <typename T, typename Alloc>
template <typename input_iterator>
stable_vector<T,Alloc>::stable_vector( input_iterator begin, input_iterator end,
                                       const allocator_type& alloc ):
v_(begin,end,alloc)
{}

template <typename T, typename Alloc>
stable_vector<T,Alloc>::stable_vector ( std::initializer_list<value_type> il,
                                         const allocator_type& alloc ):
v_(il,alloc)
{}

template <typename T, typename Alloc>
stable_vector<T,Alloc>::stable_vector( const stable_vector&  rhs,
                                      const allocator_type& alloc ):
v_(rhs,alloc)
{}

template <typename T, typename Alloc>
stable_vector<T,Alloc>::stable_vector(       stable_vector&& rhs,
                                      const allocator_type&  alloc ):
v_(rhs,alloc)
{}

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>:: begin()       noexcept ->       iterator
{ return iterator(&v_,0); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>:: begin() const noexcept -> const_iterator
{ return const_iterator(&v_,0); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::cbegin() const noexcept -> const_iterator
{ return const_iterator(&v_,0); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>:: end()       noexcept ->       iterator
{ return iterator(&v_,v_.size()); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>:: end() const noexcept -> const_iterator
{ return const_iterator(&v_,v_.size()); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::cend() const noexcept -> const_iterator
{ return const_iterator(&v_,v_.size()); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>:: rbegin()       noexcept ->       reverse_iterator
{ return reverse_iterator(end()); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>:: rbegin() const noexcept -> const_reverse_iterator
{ return const_reverse_iterator(end()); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::crbegin() const noexcept -> const_reverse_iterator
{ return const_reverse_iterator(cend()); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>:: rend()       noexcept ->       reverse_iterator
{ return reverse_iterator(begin()); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>:: rend() const noexcept -> const_reverse_iterator
{ return const_reverse_iterator(begin()); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::crend() const noexcept -> const_reverse_iterator
{ return const_reverse_iterator(cbegin()); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::empty() const noexcept -> bool
{ return v_.empty(); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::size() const noexcept -> size_type
{ return v_.size(); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::max_size() const noexcept -> size_type
{ return v_.max_size(); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::capacity() const noexcept -> size_type
{ return v_.capacity(); }

template <typename T, typename Alloc> inline
void stable_vector<T,Alloc>::reserve( size_type n )
{ v_.reserve(n); }

template <typename T, typename Alloc> inline
void stable_vector<T,Alloc>::shrink_to_fit()
{ v_.shrink_to_fit(); }

template <typename T, typename Alloc> inline
void stable_vector<T,Alloc>::resize( size_type n )
{ v_.resize(n); }

template <typename T, typename Alloc> inline
void stable_vector<T,Alloc>::resize( size_type n, const value_type& v )
{ v_.resize(n); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::operator[]( size_type n ) -> reference
{ return v_[n]; }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::operator[]( size_type n ) const -> const_reference
{ return v_[n]; }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::at( size_type n ) -> reference
{ return v_.at(n); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::at( size_type n ) const -> const_reference
{ return v_.at(n); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::front() -> reference
{ return v_.front(); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::front() const -> const_reference
{ return v_.front(); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::back() -> reference
{ return v_.back(); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::back() const -> const_reference
{ return v_.back(); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::data() noexcept -> pointer
{ return v_.data(); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::data() const noexcept -> const_pointer
{ return v_.data(); }

template <typename T, typename Alloc>
template <class input_iterator> inline
void stable_vector<T,Alloc>::assign( input_iterator first, input_iterator last )
{ v_.assign( first, last ); }

template <typename T, typename Alloc> inline
void stable_vector<T,Alloc>::assign( size_type n, const value_type& val )
{ v_.assign( n, val ); }

template <typename T, typename Alloc> inline
void stable_vector<T,Alloc>::assign( std::initializer_list<value_type> il )
{ v_.assign( il ); }

template <typename T, typename Alloc> inline
void stable_vector<T,Alloc>::push_back( const value_type&  val )
{ v_.push_back( val ); }

template <typename T, typename Alloc> inline
void stable_vector<T,Alloc>::push_back(       value_type&& val )
{ v_.push_back( val ); }

template <typename T, typename Alloc> inline
void stable_vector<T,Alloc>::pop_back()
{ v_.pop_back(); }

template <typename T, typename Alloc> inline
void stable_vector<T,Alloc>::swap( stable_vector<T,Alloc> &v )
{ v_.swap( v.v_ ); }

template <typename T, typename Alloc> inline
void stable_vector<T,Alloc>::clear() noexcept
{ v_.clear(); }

template <typename T, typename Alloc> inline
auto stable_vector<T,Alloc>::get_allocator() const noexcept -> allocator_type
{ return v_.get_allocator(); }

template <typename T, typename Alloc> inline
void swap( stable_vector<T,Alloc>& x, stable_vector<T,Alloc>& y )
{
    x.swap(y);
}

}

