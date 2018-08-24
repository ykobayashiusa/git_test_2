/*
 * Copyright (C) 2016 Matthias Kirchhart
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
#ifndef FEM_QUADRULE_ITERATOR_H
#define FEM_QUADRULE_ITERATOR_H

#include "fem/multigrid.h"
#include "geometry/quadrules.h"

namespace fem
{

template <size_t order>
class quadrule_iterator
{
public:
    using size_type         = size_t;
    using difference_type   = size_t;
    using value_type        = const geometry::quad_node;
    using reference         = const geometry::quad_node&;
    using pointer           = const geometry::quad_node*;
    using iterator_category = std::input_iterator_tag;

    using quad_node_iterator = typename tetrahedral_quadrule::const_iterator;

    const_grid_iterator<order> t;
    quad_node_iterator         n, nbegin, nend;

    bool operator==( const quadrule_iterator &rhs ) const noexcept
    { return t == rhs.t && n == rhs.n; }

    bool operator!=( const quadrule_iterator &rhs ) const noexcept
    { return t != rhs.t || n != rhs.n; }

    value_type operator* () noexcept
    { return value_type { t->Chi(n->b), t->vol_elem(n->b)*(n->w)/6 }; }

    pointer    operator->() noexcept
    { return &(this->operator*()); }

    quadrule_iterator operator++()    noexcept
    { if ( ++n == nend ) { n = nbegin; ++t; } return *this; }
        

    quadrule_iterator operator++(int) noexcept
    {
        quadrule_iterator result { *this };
        this->operator++();
        return result;
    }


    static quadrule_iterator begin( const_grid_iterator<order> grid_begin, const_grid_iterator<order>         ,
                                    quad_node_iterator         rule_begin, quad_node_iterator         rule_end ) noexcept
    { return quadrule_iterator { grid_begin, rule_begin, rule_begin, rule_end }; }

    static quadrule_iterator end  ( const_grid_iterator<order>           , const_grid_iterator<order> grid_end,
                                    quad_node_iterator         rule_begin, quad_node_iterator         rule_end ) noexcept
    { return quadrule_iterator { grid_end, rule_begin, rule_begin, rule_end }; }
};

}

#endif

