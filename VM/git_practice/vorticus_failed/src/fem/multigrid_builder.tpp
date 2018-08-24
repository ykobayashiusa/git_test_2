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

namespace fem
{

template <size_t order>
tetrahedron<order>& multigrid_builder<order>::add( multigrid<order> &mg, const nodes_t &nodes )
{
    tetrahedron<order> t;
    t.m_      = &mg;
    t.nodes_  = nodes;
    t.level_  = 0;

    triangle<order> &f0 = make_face( mg, t, 1, 2, 3 );
    triangle<order> &f1 = make_face( mg, t, 0, 2, 3 );
    triangle<order> &f2 = make_face( mg, t, 0, 1, 3 );
    triangle<order> &f3 = make_face( mg, t, 0, 1, 2 );

    t.faces_[0] = f0.id(); f0.link(t.id());
    t.faces_[1] = f1.id(); f1.link(t.id());
    t.faces_[2] = f2.id(); f2.link(t.id());
    t.faces_[3] = f3.id(); f3.link(t.id());

    // Get edges from the faces.
    t.edges_[0] = f2.edges_[2]; // (p0,p1).
    t.edges_[1] = f1.edges_[2]; // (p0,p2).
    t.edges_[2] = f0.edges_[2]; // (p1,p2).
    t.edges_[3] = f1.edges_[1]; // (p0,p3).
    t.edges_[4] = f0.edges_[1]; // (p1,p3).
    t.edges_[5] = f0.edges_[0]; // (p2,p3).

    return mg.add(t);
}

template <size_t order>
triangle<order>& multigrid_builder<order>::make_face
(
  multigrid<order> &mg, const tetrahedron<order> &t,
  uchar p0, uchar p1, uchar p2
)
{
    constexpr auto lat3d = shapefcts3d::lattice<order>();
    const point   x0 { t.nodes_[p0] }; 
    const point   x1 { t.nodes_[p1] };
    const point   x2 { t.nodes_[p2] };
    const bary3d bx0 { lat3d[p0] };
    const bary3d bx1 { lat3d[p1] };
    const bary3d bx2 { lat3d[p2] };

    triangle<order> f;
    f.m_     = &mg;
    f.level_ = 0;
    f.nodes_[0] = x0;
    f.nodes_[1] = x1;
    f.nodes_[2] = x2;
   
    constexpr auto lat2d = shapefcts2d::lattice<order>();
    triangle<order>* fptr = mg.find_face(f.id()); 
    if ( fptr == nullptr )
    {
        f.edges_[0] = make_edge(mg,t,p1,p2);
        f.edges_[1] = make_edge(mg,t,p0,p2);
        f.edges_[2] = make_edge(mg,t,p0,p1);

        // Set the higher-order DOFs.
        for ( size_t k = 3; k < lat2d.size(); ++k )
        {
            const bary2d b = lat2d[k];
            f.nodes_[k] = t.Chi( b.z0*bx0 + b.z1*bx1 + b.z2*bx2 );
        }
        return mg.add(f);
    }
    else return *fptr;
}

template <size_t order>
geoid multigrid_builder<order>::make_edge
(
  multigrid<order> &mg, const tetrahedron<order> &t,
  uchar p0, uchar p1
)
{
    constexpr auto lat3d = shapefcts3d::lattice<order>();
    const point   x0 { t.nodes_[p0] }; 
    const point   x1 { t.nodes_[p1] };
    const bary3d bx0 { lat3d[p0] };
    const bary3d bx1 { lat3d[p1] };

    edge<order> e;
    e.m_     = &mg;
    e.level_ = 0;
    e.nodes_.front() = x0;
    e.nodes_.back()  = x1;
  
    constexpr auto lat1d = shapefcts1d::lattice<order>();
    if ( ! mg.has_edge(e.id()) )
    {
        // Set the higher-order DOFs.
        for ( size_t k = 1; k < order; ++k )
        {
            const real x = lat1d[k];
            e.nodes_[k] = t.Chi( (1-x)*bx0 + x*bx1 );
        }
        mg.add(e);
    }

    return e.id();
}

}

