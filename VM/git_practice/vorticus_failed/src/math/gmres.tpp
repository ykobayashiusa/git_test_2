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

#include <vector>

#include "math/vector.h"
#include "math/matrix.h"

namespace math
{

namespace gmres_impl
{

inline
void generate_plane_rotation( real dx, real dy, real &cs, real &sn ) noexcept
{
    using std::abs;
    using std::sqrt;

    if ( dy == 0.0 )
    {
        cs = 1.0;
        sn = 0.0;
    }
    else if ( abs(dy) > abs(dx) )
    {
        real tmp = dx / dy;
        sn = 1 / sqrt( 1 + tmp*tmp );
        cs = tmp * sn;
    }
    else
    {
        real tmp = dy / dx;
        cs = 1 / sqrt( 1 + tmp*tmp );
        sn = tmp * cs;
    }
}

inline
void apply_plane_rotation( real &dx, real &dy, real cs, real sn ) noexcept
{
    real tmp = cs*dx + sn*dy;
    dy = -sn*dx + cs*dy;
    dx = tmp;
}


template <typename vec>
void update( vec &x, size_t k, math::matrix &h, math::vector &s,
             std::vector<vec> &v )
{
    math::vector y = s;

    // Solve H*y = s, H is upper triangular.
    for ( size_t i = k + 1; i-- > 0; )
    {
        y(i) /= h(i,i);
        for ( size_t j = i; j-- > 0;  )
            y(j) -= h(j,i) * y(i);
    }

    for ( size_t j = 0; j <= k; j++ )
    {
        x = x + y(j)*v[j];
    }
}

}

template <typename mat, typename vec>
real gmres( mat& A, vec& b, vec &x,
            real target_residual, size_t max_iter, bool relative )
{
    using gmres_impl::update;
    using gmres_impl::generate_plane_rotation;
    using gmres_impl::apply_plane_rotation;

    math::vector  s(max_iter+1, arma::fill::zeros ),
                 cs(max_iter+1),
                 sn(max_iter+1);
    math::matrix H(max_iter+1,max_iter+1);
    vec w;
  
    real normb = norm(b);
    if ( normb == 0.0 )
        normb = 1;

    vec r = b - A*x;
    real beta  = norm(r);
    real resid = relative ? beta/normb : beta;
    //std::cout << "GMRES-Method. Initial residual: " << resid << ".\n";
     
    if ( resid < target_residual )
        return resid;

    std::vector<vec> v( max_iter + 1 );
    v[ 0 ] = r / beta;
    s( 0 ) = beta;
    
    for ( size_t i = 0; i < max_iter; i++ )
    {
        w = A*v[i];
        for ( size_t k = 0; k <= i; k++ )
        {
            H(k,i) = dot(w,v[k]);
            w -= H(k,i)*v[k];
        }
        H(i+1,i) = norm(w);
        v[i+1]   = w / H(i+1,i);

        for ( size_t k = 0; k < i; k++ )
        {
            apply_plane_rotation( H(k,i), H(k+1,i), cs(k), sn(k) );
        }
      
        generate_plane_rotation( H(i,i), H(i+1,i), cs(i), sn(i) );
        apply_plane_rotation( H(i,i), H(i+1,i), cs(i), sn(i) );
        apply_plane_rotation( s(i), s(i+1), cs(i), sn(i) );
      
        resid = relative ? std::abs(s(i+1)) / normb : std::abs(s(i+1));
        //std::cout << "GMRES-Iteration: " << i + 1 << ", "
        //          << "Residual: " << resid << ".\n";
        if ( resid < target_residual )
        {
            update( x, i, H, s, v );
            return resid;
        }
    }

    update( x, max_iter - 1, H, s, v );
    r = b - A*x;
    beta = norm(r);
    resid = relative ? beta / normb : beta;
    return resid;
}

template <typename precond, typename mat, typename vec>
real gmres( precond &Q, mat& A, vec& b, vec &x,
            real target_residual, size_t max_iter, bool relative )
{
    using gmres_impl::update;
    using gmres_impl::generate_plane_rotation;
    using gmres_impl::apply_plane_rotation;

    math::vector  s(max_iter+1, arma::fill::zeros ),
                 cs(max_iter+1),
                 sn(max_iter+1);
    math::matrix H(max_iter+1,max_iter+1);
    vec w;
  
    real normb = norm(Q*b);
    if ( normb == 0.0 )
        normb = 1;

    vec r = Q*(b - A*x);
    real beta  = norm(r);
    real resid = relative ? beta/normb : beta;
    std::cout << "GMRES-Method. Initial residual: " << resid << ".\n";
    if ( resid < target_residual )
        return resid;

    std::vector<vec> v( max_iter + 1 );
    v[ 0 ] = r / beta;
    s( 0 ) = beta;
    
    for ( size_t i = 0; i < max_iter; i++ )
    {
        w = Q*(A*v[i]);
        for ( size_t k = 0; k <= i; k++ )
        {
            H(k,i) = dot(w,v[k]);
            w -= H(k,i)*v[k];
        }
        H(i+1,i) = norm(w);
        v[i+1]   = w / H(i+1,i);

        for ( size_t k = 0; k < i; k++ )
        {
            apply_plane_rotation( H(k,i), H(k+1,i), cs(k), sn(k) );
        }
      
        generate_plane_rotation( H(i,i), H(i+1,i), cs(i), sn(i) );
        apply_plane_rotation( H(i,i), H(i+1,i), cs(i), sn(i) );
        apply_plane_rotation( s(i), s(i+1), cs(i), sn(i) );
      
        resid = relative ? std::abs(s(i+1)) / normb : std::abs(s(i+1));
        std::cout << "GMRES-Iteration: " << i + 1 << ", "
                  << "Residual: " << resid << ".\n";
        if ( resid < target_residual )
        {
            update( x, i, H, s, v );
            return resid;
        }
    }

    update( x, max_iter - 1, H, s, v );
    r = Q*(b - A*x);
    beta = norm(r);
    resid = relative ? beta / normb : beta;
    return resid;
}

}

