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
#ifndef FEM_MULTIGRID_BUILDER_H
#define FEM_MULTIGRID_BUILDER_H

#include "fem/multigrid.h"

namespace fem
{

template <size_t order>
class multigrid_builder
{
public:
    virtual ~multigrid_builder() = default;

protected:
    using nodes_t = shapefcts3d::coeffs<order,point>;
    tetrahedron<order>& add( multigrid<order> &mg, const nodes_t &nodes );

private:
    triangle<order>& make_face( multigrid<order> &mg, const tetrahedron<order> &t,
                                uchar p0, uchar p1, uchar p2 );
    geoid            make_edge( multigrid<order> &mg, const tetrahedron<order> &t,
                                uchar p0, uchar p1 );
};

}


#include "fem/multigrid_builder.tpp"
#endif

