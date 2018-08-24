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

namespace math
{

class hash_matrix::handle
{
private:
    using index = hash_matrix::index;
    using map   = decltype(hash_matrix::values_);

    index  pos;
    map  &data;

public:
    handle( index p, map &d ): pos { p }, data { d } {}

    handle() = delete;
    handle( const handle&  rhs ) = delete;
    handle(       handle&& rhs ) = default;
    handle& operator=( const handle&  rhs ) = delete;
    handle& operator=(       handle&& rhs ) = delete;

    operator real() const noexcept
    { auto it = data.find(pos); return (it!=data.end()) ? it->second : 0; }

    real& operator =( real rhs ) { data[pos]  = rhs; return data[pos]; }
    void  operator+=( real rhs ) { data[pos] += rhs; }
    void  operator-=( real rhs ) { data[pos] -= rhs; }
    void  operator*=( real rhs ) noexcept
    { auto it = data.find(pos); if (it!=data.end()) it->second *= rhs; } 
    void  operator/=( real rhs ) noexcept
    { auto it = data.find(pos); if (it!=data.end()) it->second /= rhs; } 
};

inline
hash_matrix::hash_matrix( size_t n_rows, size_t n_cols ):
rows_ { n_rows }, cols_ { n_cols }
{}

inline
hash_matrix::handle hash_matrix::operator()( size_t row, size_t col ) noexcept
{
    assert( row < rows_ && col < cols_ ); 
    return handle { index { row, col }, values_ };
}

inline
real hash_matrix::operator()( size_t row, size_t col ) const noexcept
{
    assert( row < rows_ && col < cols_ ); 
    auto entry = values_.find( index { row, col } );
    return ( entry == values_.end() ) ? 0 : entry->second;
}

inline
void hash_matrix::set_value( size_t row, size_t col, real value )
{
    (*this)(row,col) = value;
}

inline
void hash_matrix::del_value( size_t row, size_t col )
{
    assert( row < rows_ && col < cols_ ); 
    values_.erase( index { row, col } );
}

inline
void hash_matrix::operator*=( real rhs ) noexcept
{
    for ( auto &entry: values_ )
        entry.second *= rhs;
}

inline
void hash_matrix::operator/=( real rhs ) noexcept
{
    for ( auto &entry: values_ )
        entry.second /= rhs;
}

inline
void hash_matrix::operator+=( const hash_matrix &rhs )
{
    assert( rows_ == rhs.rows_ && cols_ == rhs.cols_ );
    for ( auto &entry: rhs.values_ )
        values_[ entry.first ] += entry.second;
}

inline
void hash_matrix::operator-=( const hash_matrix &rhs )
{
    assert( rows_ == rhs.rows_ && cols_ == rhs.cols_ );
    for ( auto &entry: rhs.values_ )
        values_[ entry.first ] -= entry.second;
}

inline
void hash_matrix::clear() noexcept
{
    values_.clear();
}

inline
size_t hash_matrix::rows() const noexcept
{
    return rows_;
}

inline
size_t hash_matrix::cols() const noexcept
{
    return cols_;
}

inline
size_t hash_matrix::nnz() const noexcept
{
    return values_.size();
}

inline
void hash_matrix::resize( size_t n_rows, size_t n_cols )
{
    const bool enlarged { n_rows >= rows_ && n_cols >= cols_ };

    rows_ = n_rows;
    cols_ = n_cols;

    if ( enlarged )
        return;

    auto it = values_.cbegin();
    while ( it != values_.end() )
    {
        index idx = it->first;
        if ( idx.row < rows_ && idx.col < cols_ )
        {
            ++it;
        }
        else
        {
            values_.erase( it );
        }
    }
}

inline
real hash_matrix::infty_norm() const
{
    std::vector<real> row_values( rows_, 0 );
    for ( const auto& entry: values_ )
    {
        size_t row = entry.first.row;
        row_values[ row ] += std::abs( entry.second );
    }
    return *std::max_element( row_values.begin(), row_values.end() );
}

inline
auto hash_matrix::begin() noexcept -> iterator
{
    return values_.begin();
}

inline
auto hash_matrix::end() noexcept -> iterator
{
    return values_.end();
}

inline
auto hash_matrix::begin() const noexcept -> const_iterator
{
    return values_.begin();
}

inline
auto hash_matrix::end() const noexcept -> const_iterator
{
    return values_.end();
}

inline
auto hash_matrix::cbegin() const noexcept -> const_iterator
{
    return values_.cbegin();
}

inline
auto hash_matrix::cend() const noexcept -> const_iterator
{
    return values_.cend();
}

inline
vector hash_matrix::operator*( const vector& rhs ) const
{
    assert( cols_ == rhs.size() );
    vector result( rows_, arma::fill::zeros );
    for ( const value_type& entry: values_ )
    {
        const size_t row = entry.first.row;
        const size_t col = entry.first.col;
        const real   val = entry.second;

        result(row) += val*rhs(col);
    }
    return std::move(result);
}



struct crs_matrix::handle
{
    crs_matrix &A;
    size_t row, col;

    operator real() const noexcept { return ((const crs_matrix&) A)(row,col); }

    real operator*( real rhs ) const noexcept { return A(row,col)*rhs; }
    real operator/( real rhs ) const noexcept { return A(row,col)/rhs; }
    real operator+( real rhs ) const noexcept { return A(row,col)+rhs; }
    real operator-( real rhs ) const noexcept { return A(row,col)-rhs; }

    real operator=( real rhs ) { A.set_value(row,col,rhs); return rhs; }

    void operator*=( real rhs ) noexcept
    {
        value_iterator it = A.find_value(row,col);
        if ( it != A.values_.end() )
            (*it) *= rhs;
    }

    void operator/=( real rhs ) noexcept
    {
        value_iterator it = A.find_value(row,col);
        if ( it != A.values_.end() )
            (*it) /= rhs;
    }

    void operator+=( real rhs )
    {
        value_iterator it = A.find_value(row,col);
        if ( it != A.values_.end() ) (*it) += rhs;
        else A.set_value( row, col, rhs );
    }

    void operator-=( real rhs )
    {
        value_iterator it = A.find_value(row,col);
        if ( it != A.values_.end() ) (*it) -= rhs;
        else A.set_value( row, col, -rhs );
    }
};

inline
crs_matrix::crs_matrix( size_t p_rows, size_t p_cols ):
rows_ { p_rows }, cols_ { p_cols },
row_ptr_( rows_+ 1 ), col_idx_(0), values_(0)
{}

inline
crs_matrix::crs_matrix( const hash_matrix& A ):
rows_ { A.rows() }, cols_ { A.cols() },
row_ptr_( rows_ + 1 ), col_idx_( A.nnz() ), values_( A.nnz() )
{
    using entry = std::pair< hash_matrix::index, real >;
    auto cmp = []( entry lhs, entry rhs ) noexcept -> bool
    {
        if ( lhs.first.row == rhs.first.row ) return lhs.first.col < rhs.first.col;
        else                                  return lhs.first.row < rhs.first.row;
    };

    std::vector<entry> entries( A.begin(), A.end() );
    std::sort( entries.begin(), entries.end(), cmp );

    size_t nnz_count = 0;
    for ( size_t row = 0; row <= rows_; ++row )
    {
        row_ptr_[ row ] = nnz_count;
        while ( nnz_count < A.nnz() && entries[ nnz_count ].first.row == row )
        {
            col_idx_[ nnz_count ] = entries[ nnz_count ].first.col;
            values_ [ nnz_count ] = entries[ nnz_count ].second;
            ++nnz_count;
        }
    } 
}

template <typename matrix_entry_iterator>
crs_matrix::crs_matrix( size_t p_rows, size_t p_cols,
                        matrix_entry_iterator begin, matrix_entry_iterator end ):
rows_ { p_rows }, cols_ { p_cols },
row_ptr_( rows_ + 1, 0 ),
col_idx_( end - begin ),
values_ ( end - begin )
{
          std::sort  (begin,end,row_major_cmp);
    end = std::unique(begin,end,row_major_cmp);

    size_t nnz_count = 0;
    for ( size_t row = 0; row <= rows_; ++row )
    {
        assert( begin->row < rows_ );
        row_ptr_[ row ] = nnz_count;
        while ( begin != end && begin->row == row )
        {
            assert( begin->col < cols_ );
            col_idx_[ nnz_count ] = begin->col;
             values_[ nnz_count ] = begin->value;
            ++nnz_count; ++begin;
        }
    } 
}

inline
size_t crs_matrix::rows() const noexcept
{
    return rows_;
}

inline
size_t crs_matrix::cols() const noexcept
{
    return cols_;
}

inline
size_t crs_matrix::nnz() const noexcept
{
    return values_.size();
}

inline
crs_matrix::handle crs_matrix::operator()( size_t row, size_t col )
{
    return handle { *this, row, col };
}

inline
real crs_matrix::operator()( size_t row, size_t col ) const noexcept
{
    using iter = std::vector<size_t>::const_iterator;

    assert( row < rows_ && col < cols_ );
    iter begin = col_idx_.cbegin() + row_ptr_[ row ];
    iter end   = col_idx_.cbegin() + row_ptr_[ row + 1 ];

    iter pos = std::lower_bound( begin, end, col );

    if ( pos != end && *pos == col ) return values_[ pos - col_idx_.begin() ];
    else return 0;
}

inline
void crs_matrix::set_value( size_t row, size_t col, real value )
{
    using iter = std::vector<size_t>::iterator;
    assert( row < rows_ && col < cols_ );

    iter begin = col_idx_.begin() + row_ptr_[ row ];
    iter end   = col_idx_.begin() + row_ptr_[ row + 1 ];

    iter pos = std::lower_bound( begin, end, col );

    if ( pos != end && *pos == col )
    {
        values_[ pos - col_idx_.begin() ] = value;
    }
    else
    {
        auto value_pos = values_.begin() + ( pos - col_idx_.begin() );
        col_idx_.insert( pos, col );
         values_.insert( value_pos, value );
        for ( size_t i = row + 1; i <= rows_; ++i )
            ++row_ptr_[ i ];
    }
}

inline
void crs_matrix::del_value( size_t row, size_t col ) noexcept
{
    using iter = std::vector<size_t>::iterator;

    assert( row < rows_ && col < cols_ );
    iter begin = col_idx_.begin() + row_ptr_[ row ];
    iter end   = col_idx_.begin() + row_ptr_[ row + 1 ];
    iter pos = std::lower_bound( begin, end, col );

    if ( pos != end && *pos == col )
    {
         values_.erase( values_.begin() + (pos - col_idx_.begin()));
        col_idx_.erase( pos );
        for ( size_t i = row + 1;  i <= rows_; ++i )
            --row_ptr_[ i ];
    }
}

inline
void crs_matrix::del_row( size_t row ) noexcept
{
    assert( row < rows_ );
    auto  begin = col_idx_.begin() + row_ptr_[ row ];
    auto vbegin =  values_.begin() + row_ptr_[ row ];
    auto  end   = col_idx_.begin() + row_ptr_[ row + 1 ];
    auto vend   =  values_.begin() + row_ptr_[ row + 1 ];
    size_t n_elems = end - begin;

    col_idx_.erase(  begin,  end );
     values_.erase( vbegin, vend );
    for ( size_t i = row + 1; i <= rows_; ++i )
        row_ptr_[ i ] -= n_elems;
}

inline
void crs_matrix::clear() noexcept
{
    row_ptr_.assign( rows_ + 1, 0 );
    col_idx_.clear();
     values_.clear();
}

inline
void crs_matrix::compress()
{
    for ( size_t row = 0; row < rows_; ++row )
    {
        for ( size_t i = row_ptr_[row]; i < row_ptr_[row+1]; ++i )
        {
            if ( values_[i] == 0 )
            {
                 values_.erase(  values_.begin() + i );
                col_idx_.erase( col_idx_.begin() + i );
                for ( size_t j = row + 1;  j <= rows_; ++j )
                    --row_ptr_[ j ];
            }
        }
    }

    row_ptr_.shrink_to_fit();
    col_idx_.shrink_to_fit();
     values_.shrink_to_fit();
}

inline
vector crs_matrix::operator*( const vector& rhs ) const
{
    assert( cols_ == rhs.size() );
    vector result( rows_, arma::fill::zeros );

    #pragma omp parallel for schedule(dynamic)
    for ( size_t i = 0; i < rows_; ++i )
    {
        double tmp = 0;
        for ( size_t j = row_ptr_[ i ]; j < row_ptr_[ i + 1 ]; ++j )
        {
            tmp += values_[ j ]*rhs( col_idx_[ j ] );
        }
        result(i) += tmp;
    }
    return result;
}

inline
void crs_matrix::operator*=( real rhs ) noexcept
{
    #pragma omp parallel for
    for ( size_t i = 0; i < values_.size(); ++i )
        values_[ i ] *= rhs;
}

inline
void crs_matrix::operator/=( real rhs ) noexcept
{
    #pragma omp parallel for
    for ( size_t i = 0; i < values_.size(); ++i )
        values_[ i ] /= rhs;
}

inline
crs_matrix::value_iterator crs_matrix::find_value( size_t row, size_t col ) noexcept
{
    assert( row < rows_ && col < cols_ );
    auto begin = col_idx_.begin() + row_ptr_[ row ];
    auto end   = col_idx_.begin() + row_ptr_[ row + 1 ];
    auto pos = std::lower_bound( begin, end, col );

    if ( pos != end && *pos == col ) return values_.begin() + (pos-begin);
    else return values_.end();
}

inline
crs_matrix::const_value_iterator crs_matrix::find_value( size_t row, size_t col ) const noexcept
{
    assert( row < rows_ && col < cols_ );
    auto begin = col_idx_.cbegin() + row_ptr_[ row ];
    auto end   = col_idx_.cbegin() + row_ptr_[ row + 1 ];
    auto pos = std::lower_bound( begin, end, col );

    if ( pos != end && *pos == col ) return values_.cbegin() + (pos-begin);
    else return values_.cend();
}

}

