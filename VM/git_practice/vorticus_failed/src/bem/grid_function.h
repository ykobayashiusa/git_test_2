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
 * details.  *
 * You should have received a copy of the GNU General Public License along with
 * vorticus; see the file COPYING.  If not see http://www.gnu.org/licenses.
 */
#ifndef BEM_GRID_FUNCTION_H
#define BEM_GRID_FUNCTION_H

#include "bem/multigrid.h"

#include <boost/functional/hash.hpp>

namespace bem
{

template <typename T, uint gorder>
class grid_function_base
{
public:
    virtual T eval( point pos ) const = 0;
    virtual T eval( const triangle<gorder>& t, bary pos ) const = 0;
};

template <uint gorder, uint aorder>
class prolongation;

template <typename T, uint gorder, uint aorder>
class grid_function: public grid_function_base<T,gorder>
{
public:
    friend class prolongation<gorder,aorder>;

    template <uint o>
    using grid = typename multigrid<o>::grid;

    using       coords_t = std::array<point,lattice::size<aorder>()>;
    using       iterator = typename std::unordered_map<point,T>::      iterator;
    using const_iterator = typename std::unordered_map<point,T>::const_iterator;
    using      size_type = typename std::unordered_map<point,T>::size_type;

    grid_function( const grid<gorder>& g );
    grid_function( const grid_function&  rhs ) = default;
    grid_function(       grid_function&& rhs ) = default;
    grid_function& operator=( const grid_function&  rhs ) = default;
    grid_function& operator=(       grid_function&& rhs ) = default;

          T& operator()( point pos ) noexcept;
    const T& operator()( point pos ) const noexcept;

    T eval( point pos ) const override;
    T eval( const triangle<gorder>& t, bary pos ) const override;

    iterator begin() noexcept;
    iterator   end() noexcept;

    const_iterator begin() const noexcept;
    const_iterator   end() const noexcept;
    
    size_type size() const noexcept;

    unsigned char       get_level() const noexcept { return level_; }
    coords_t            get_dofs( const triangle<gorder>& t ) const noexcept;
    

    void operator+=( const grid_function& rhs );
    void operator-=( const grid_function& rhs );
    void operator*=( const real rhs );
    void operator/=( const real rhs );
    grid_function& operator=( const T& rhs );


private:
    unsigned char level_;
    std::unordered_map<point,T>  values_;
    std::unordered_map<geoid,coords_t> coords_;
};

template <uint gorder, uint aorder>
class prolongation
{
public:
    prolongation() = default;
    prolongation( const typename multigrid<gorder>::grid& coarse_grid,
                  const typename multigrid<gorder>::grid&   fine_grid );

    template <typename T>
    void prolongate( const grid_function<T,gorder,aorder>& coarse_function,
                           grid_function<T,gorder,aorder>&   fine_function ) const;

    template <typename T>
    void   restrict(       grid_function<T,gorder,aorder>& coarse_function,
                     const grid_function<T,gorder,aorder>&   fine_function ) const;
    
private:
    using idx_t = std::pair<point,point>;
    struct idx_hash
    {
        size_t operator()( idx_t pp ) const noexcept
        {
            using boost::hash_combine;

            std::hash<point> hasher;
            size_t result { 0 };
            hash_combine( result, hasher(pp.first)  );
            hash_combine( result, hasher(pp.second) );
            return result;
        }
    };

    using scalar = grid_function<real,gorder,aorder>;
    void treat_triangle( const triangle<gorder>& t, const scalar &coarse, const scalar &fine );
    void treat_no_ref          ( const triangle<gorder>& t, const scalar &coarse, const scalar &fine );
    void treat_single_irreg_ref( const triangle<gorder>& t, const scalar &coarse, const scalar &fine );
    void treat_double_irreg_ref( const triangle<gorder>& t, const scalar &coarse, const scalar &fine );
    void treat_regular_ref     ( const triangle<gorder>& t, const scalar &coarse, const scalar &fine );

    std::unordered_map<idx_t,real,idx_hash> matrix_ {};
};

}

#include "bem/grid_function.tpp"
#endif

