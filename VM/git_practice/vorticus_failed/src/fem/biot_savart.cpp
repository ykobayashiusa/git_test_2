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
#include "fem/biot_savart.h"

#include "math/constants.h" 
#include "geometry/point.h"

#include <iostream>

namespace fem
{

using geometry::point;

point biot_savart( const tetrahedron<1> &t, const shapefcts3d::coeffs<1,point> &coeffs, const point P )
{
    using shapefcts3d:: N;
    using shapefcts3d::dNdxi;
    using shapefcts3d::dNdeta;
    using shapefcts3d::dNdzeta;
    using shapefcts3d::apply;
    using math::eps;
    using std::min;

    // Compute the gradient and curl of omega.
    const point dxi   = apply<1>( coeffs, dNdxi  <1>( bary3d {} ) );
    const point deta  = apply<1>( coeffs, dNdeta <1>( bary3d {} ) );
    const point dzeta = apply<1>( coeffs, dNdzeta<1>( bary3d {} ) );
    const tensor grad_omega = ( tensor { dxi.x, deta.x, dzeta.x,
                                         dxi.y, deta.y, dzeta.y,
                                         dxi.z, deta.z, dzeta.z } ) *
                              t.dChi( bary3d {} ).inv();


    const point rot_omega { grad_omega(2,1) - grad_omega(1,2),
                            grad_omega(0,2) - grad_omega(2,0),
                            grad_omega(1,0) - grad_omega(0,1) };
   
    const bary3d bar { t.ChiAffInv(P) }; 
    const real min_bar { min( min(bar.z0,bar.z1), min(bar.z2,bar.z3) ) };

    real alpha;
         if ( min_bar >  0 ) alpha = 4*math::pi;
    else if ( min_bar == 0 ) alpha = 2*math::pi;
    else                     alpha = 0;
    
    // Corner points making up the individual faces. Last point
    // equals first point, so that edge vectors can be computed
    // easily.
    const point x0 = t.get_node(0); const point x1 = t.get_node(1);
    const point x2 = t.get_node(2); const point x3 = t.get_node(3);
    const point  x[4][4]
    {
        { x1, x2, x3, x1 },
        { x0, x3, x2, x0 },
        { x0, x1, x3, x0 },
        { x0, x2, x1, x0 }
    };

    // The following two arrays are to avoid the unnecessary recomputation
    // of distances for the edges, as they are independent of the direction
    // the edges are integrated over.
    // Distances from corners to the field point.
    const real r0 = (P-x0).r(); const real r1 = (P-x1).r();
    const real r2 = (P-x2).r(); const real r3 = (P-x3).r();
    const real r[4][4]
    {
        { r1, r2, r3, r1 },
        { r0, r3, r2, r0 },
        { r0, r1, r3, r0 }, 
        { r0, r2, r1, r0 }
    };

    // Lengths of the edges.
    const real l0 = (x1-x0).r(); const real l1 = (x2-x0).r();
    const real l2 = (x2-x1).r(); const real l3 = (x3-x0).r();
    const real l4 = (x3-x1).r(); const real l5 = (x3-x2).r();
    const real l[4][3]
    {
        { l2, l5, l4 },
        { l3, l5, l1 },
        { l0, l4, l3 },
        { l1, l2, l0 }
    };

    // Log terms.
    real E0 = ( (r1+l0)*(r1+l0) - r0*r0) + eps; E0 /= (-(r0-l0)*(r0-l0) + r1*r1) + eps; E0 = std::log(E0);
    real E1 = ( (r2+l1)*(r2+l1) - r0*r0) + eps; E1 /= (-(r0-l1)*(r0-l1) + r2*r2) + eps; E1 = std::log(E1);
    real E2 = ( (r2+l2)*(r2+l2) - r1*r1) + eps; E2 /= (-(r1-l2)*(r1-l2) + r2*r2) + eps; E2 = std::log(E2);
    real E3 = ( (r3+l3)*(r3+l3) - r0*r0) + eps; E3 /= (-(r0-l3)*(r0-l3) + r3*r3) + eps; E3 = std::log(E3);
    real E4 = ( (r3+l4)*(r3+l4) - r1*r1) + eps; E4 /= (-(r1-l4)*(r1-l4) + r3*r3) + eps; E4 = std::log(E4);
    real E5 = ( (r3+l5)*(r3+l5) - r2*r2) + eps; E5 /= (-(r2-l5)*(r2-l5) + r3*r3) + eps; E5 = std::log(E5);

    const real E[4][3]
    {
        { E2, E5, E4 },
        { E3, E5, E1 },
        { E0, E4, E3 },
        { E1, E2, E0 }
    };

    const point tcentre { t.centre() };

    point result {};
    for ( size_t j = 0; j < 4; ++j )
    {
        point n     { cross_prod( (x[j][1]-x[j][0]), (x[j][2]-x[j][0]) ) }; n    /= n   .r();
        point e_xi  { n.y*n.y + n.z*n.z, -n.x*n.y, -n.x*n.z };              e_xi /= e_xi.r();
        point e_eta { cross_prod( n, e_xi ) }; 

        point R0 = P - x[j][0];
        real xr  = scal_prod(R0,e_xi);
        real zr  = scal_prod(R0,n);
        real a0  = x[j][0].x + xr*e_xi.x;
        real a1  = e_xi.x;

        // Specially handle the case when n points in the x-direction.
        if ( n.y == 0 && n.z == 0 )
        {
            e_xi = x[j][1] - x[j][0];      e_xi  /= e_xi .r();
            e_eta = cross_prod( n, e_xi ); e_eta /= e_eta.r();
            xr = scal_prod(R0,e_xi);
            zr = scal_prod(R0,n);
            a0 = x[j][0].x;
            a1 = 0;
        }

        int sign = 1;
        real e   = zr;
        if ( e < 0 ) { e = -e; sign = -1; }

        point omega = ( j == 0 ) ? coeffs[1] : coeffs[0];
        point t0 = cross_prod( omega, n );
        point t1 = cross_prod( grad_omega*e_xi,  n );
        point t2 = cross_prod( grad_omega*e_eta, n );
        t0 += t1*scal_prod( e_xi,  R0 ) +
              t2*scal_prod( e_eta, R0 );

        point L;
        real phi_sigma0 = 0; real phi_sigma1 = 0;
        real phi_mu0    = 0; real phi_mu1    = 0; real phi_mu2   = 0;
        for ( size_t i = 0; i < 3; ++i )
        {
            const real li = l[j][i];
            const point s = (x[j][i+1] - x[j][i]) / li;
            
            const point Ri  = x[j][i] - P;
            const real  ri  = r[j][i];  
            const real  ri1 = r[j][i+1];
            const real  xd  = -scal_prod( Ri, s );
            const real  yd2 = ( Ri + xd*s ).r2() + eps;
            const real  b   = scal_prod( cross_prod(n,Ri), s );
            const real sxi  = scal_prod(s,e_xi);
            const real seta = scal_prod(s,e_eta);
            const real xi_i = scal_prod( x[j][i] - x[j][0], e_xi );

            const real sq = sqrt( yd2 - e*e ); 
            real K1_2 {};
            
            const real H  = sq*( yd2*li + e*(li-xd)*ri + e*xd*ri1 ) / ( yd2*(ri+e)*(ri1+e) );
            const real F1 = (yd2 + e*ri )/(e+ri );
            const real F2 = (yd2 + e*ri1)/(e+ri1);
            const real F  = (F1*F1 + F2*F2);
            real beta;
            if ( (xd*(li-xd)) <= 0 || F > yd2 ) beta =             asin(H);
            else if ( xd > 0 )                  beta =  math::pi - asin(H);
            else                                beta = -math::pi - asin(H);
            K1_2 = (e/sq)*beta;
            if ( ! std::isnormal(K1_2) ) K1_2 = 0;

            const real Ei = E[j][i];
            const real K1 = Ei - K1_2;
            const real K2 = .5*((li-xd)*ri1 + xd*ri + yd2*Ei);
            const real K3  = (xi_i - xr)*Ei + sxi*(ri1-ri+xd*Ei);

            L += b*t0*K1 + (seta*t1 - sxi*t2)*K2;
            phi_sigma0 += b*K1;
            phi_sigma1 += seta*K2;
                 if ( sign ==  1 && e != 0 ) phi_mu0 += b*((Ei-K1)/e);
            else if ( sign == -1 && e != 0 ) phi_mu0 -= b*((Ei-K1)/e);
            phi_mu1 += -zr*seta*Ei;
            phi_mu2 += -seta*K3;
        }
        phi_mu2 = zr*( phi_sigma0 + phi_mu2 );

        real Ksigma, Kmu;
        Ksigma =  n.x*( a0*phi_sigma0 + a1*phi_sigma1 );

        if ( a1 != 0 ) Kmu = .5*a0*a0*phi_mu0 + a0*a1*phi_mu1 + .5*a1*a1*phi_mu2;
        else           Kmu = .5*a0*a0*phi_mu0;

        real K = Ksigma - Kmu;

        if ( scal_prod( x[j][0] - tcentre, n ) < 0 ) { result -= L + K*rot_omega; }
        else                                         { result += L + K*rot_omega; }
    }
    return (result - rot_omega*alpha*P.x*P.x/2) / (4*math::pi);
}

}

