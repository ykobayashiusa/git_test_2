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
#ifndef MATH_SPMAT_H
#define MATH_SPMAT_H

#include "types.h"
#include "math/vector.h"

#include <cassert>
#include <armadillo>
#include <unordered_map>
#include <boost/functional/hash.hpp>

namespace math
{

class hash_matrix
{
public:
    struct index
    {
        size_t row, col;
        constexpr bool operator==( index rhs ) const noexcept
        {
            return (row == rhs.row) && (col == rhs.col); 
        }
    };

private:
    class handle;

    struct index_hash
    {
        size_t operator()( index i ) const noexcept
        {
            using boost::hash_combine;

            std::hash<size_t> hasher;
            size_t result { 0 };
            hash_combine( result, hasher(i.row) );
            hash_combine( result, hasher(i.col) );
            return result;
        }
    };

    size_t rows_ {0}, cols_ {0};
    std::unordered_map<index,real,index_hash> values_;


public:
    using size_type = size_t;
    using const_iterator = typename decltype( values_ )::const_iterator;
    using       iterator = typename decltype( values_ )::      iterator;
    using value_type     = typename decltype( values_ )::value_type;

    hash_matrix() = default;
    hash_matrix( const hash_matrix&  ) = default;
    hash_matrix(       hash_matrix&& ) = default;
    hash_matrix( size_t n_rows, size_t n_cols );
    hash_matrix& operator=( const hash_matrix&  ) = default;
    hash_matrix& operator=(       hash_matrix&& ) = default;

    handle  operator()( size_t row, size_t col )       noexcept;
    real    operator()( size_t row, size_t col ) const noexcept;
    void  set_value ( size_t row, size_t col, real value );
    void  del_value ( size_t row, size_t col );

    void operator*=( real rhs ) noexcept;
    void operator/=( real rhs ) noexcept;

    void operator+=( const hash_matrix &rhs );
    void operator-=( const hash_matrix &rhs );

    void clear() noexcept;

    size_t rows() const noexcept;
    size_t cols() const noexcept;
    size_t nnz()  const noexcept;

    void resize( size_t n_rows, size_t n_cols );

    real infty_norm() const;

          iterator  begin()       noexcept;
          iterator    end()       noexcept;
    const_iterator  begin() const noexcept;
    const_iterator    end() const noexcept;
    const_iterator cbegin() const noexcept;
    const_iterator   cend() const noexcept;

    vector operator*( const vector &rhs ) const; 
};



struct matrix_entry
{
    size_t row, col;
    real   value;
};

constexpr
bool row_major_cmp( matrix_entry lhs, matrix_entry rhs )
{
    return ( lhs.row != rhs.row ) ? lhs.row < rhs.row : lhs.col < rhs.col;
}

constexpr
bool col_major_cmp( matrix_entry lhs, matrix_entry rhs )
{
    return ( lhs.col != rhs.col ) ? lhs.col < rhs.col : lhs.row < rhs.row;
}


struct crs_matrix
{
public:
    struct handle;
    using       value_iterator = typename std::vector<real>  ::      iterator;
    using const_value_iterator = typename std::vector<real>  ::const_iterator;
    using         col_iterator = typename std::vector<size_t>::      iterator;
    using   const_col_iterator = typename std::vector<size_t>::const_iterator;
    using         row_iterator = typename std::vector<size_t>::      iterator;
    using   const_row_iterator = typename std::vector<size_t>::const_iterator;

    crs_matrix() = default;
    crs_matrix( const crs_matrix&  ) = default;
    crs_matrix(       crs_matrix&& ) = default;
    crs_matrix& operator=( const crs_matrix&  ) = default;
    crs_matrix& operator=(       crs_matrix&& ) = default;
   
    crs_matrix( const hash_matrix& A );
    crs_matrix( size_t p_rows, size_t p_cols );

    template <typename matrix_entry_iterator>
    crs_matrix( size_t p_rows, size_t p_cols,
                matrix_entry_iterator begin, matrix_entry_iterator end );

    size_t rows() const noexcept;
    size_t cols() const noexcept;
    size_t nnz()  const noexcept;

    handle operator()( size_t row, size_t col );
    real   operator()( size_t row, size_t col ) const noexcept;
    void  set_value( size_t row, size_t col, real value );
    void  del_value( size_t row, size_t col ) noexcept;
    void  del_row  ( size_t row ) noexcept;
    void  clear    () noexcept;

    void compress();

    vector operator*( const vector &rhs ) const; 

    void operator*=( real rhs ) noexcept;
    void operator/=( real rhs ) noexcept;


          value_iterator find_value( size_t row, size_t col )       noexcept;
    const_value_iterator find_value( size_t row, size_t col ) const noexcept;

          value_iterator values_begin ()       noexcept { return values_.begin(); }
          value_iterator values_end   ()       noexcept { return values_.end  (); }
    const_value_iterator values_begin () const noexcept { return values_.begin(); }
    const_value_iterator values_end   () const noexcept { return values_.end  (); }
    const_value_iterator values_cbegin() const noexcept { return values_.begin(); }
    const_value_iterator values_cend  () const noexcept { return values_.end  (); }

    const_col_iterator colidx_begin () const noexcept { return col_idx_.begin(); }
    const_col_iterator colidx_end   () const noexcept { return col_idx_.end  (); }
    const_col_iterator colidx_cbegin() const noexcept { return col_idx_.begin(); }
    const_col_iterator colidx_cend  () const noexcept { return col_idx_.end  (); }
    
    const_row_iterator rowptr_begin () const noexcept { return row_ptr_.begin(); }
    const_row_iterator rowptr_end   () const noexcept { return row_ptr_.end  (); }
    const_row_iterator rowptr_cbegin() const noexcept { return row_ptr_.begin(); }
    const_row_iterator rowptr_cend  () const noexcept { return row_ptr_.end  (); }

//private:
    size_t rows_ {0}, cols_ {0};
    std::vector<size_t> row_ptr_;
    std::vector<size_t> col_idx_;
    std::vector<real>   values_;
};

}

#include "math/spmat.tpp"
#endif

