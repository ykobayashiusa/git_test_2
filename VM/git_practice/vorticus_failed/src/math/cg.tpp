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

template <typename mat, typename vec>
real cg( mat& A, vec& b, vec &x,
         real target_residual, size_t max_iter, bool relative )
{
    vec r = b - A*x;
    vec p = r;
    vec Ap = A*p;

    real norm_r = norm(r);
    real norm_b = norm(b);
    real resid  = relative ? norm_r/norm_b : norm_r;

    if ( norm_b == 0 )
    {
        x.fill(0);
        resid = 0;
    }

    std::cout << "CG: Initial residual: " << resid << ".\n";
    std::cout.flush();
    if ( resid < target_residual )
        return resid;
    for ( size_t k = 0; k < max_iter; ++k )
    {
        real rr    = dot(r,r);
        real alpha = rr/dot(p,Ap);
        x = x + alpha*p;
        r = r - alpha*Ap;

        norm_r = norm(r);
        resid  = relative ? norm_r/norm_b : norm_r;
        if ( (k + 1) % 10 == 0 )
        {
            std::cout << "CG: Iteration:  " << std::setw(4) << k + 1 << ", "
                      << "Residual:  " <<  std::setw(12) << std::scientific << resid << ".\n";
            std::cout.flush();
        }
        if ( resid <= target_residual )
        {
            std::cout << "CG: Target residual "  << target_residual
                      << " reached after " << k + 1 << " iterations.\n";
            std::cout.flush();
            return resid;
        }

        real beta = dot(r,r)/rr;
        p = r + beta*p;
        Ap = A*p;
    }

    return resid;
}

template <typename precond, typename mat, typename vec>
real cg( precond &Q, mat& A, vec& b, vec &x,
         real target_residual, size_t max_iter, bool relative )
{
    vec r  = b - A*x;
    vec z  = Q*r;
    vec p  = z;
    vec Ap;

    real norm_r = norm(r);
    real norm_b = norm(b);
    real resid  = relative ? norm_r/norm_b : norm_r;

    if ( norm_b == 0 )
    {
        x.fill(0);
        resid = 0;
    }

    std::cout << "CG: Initial residual: " << resid << ".\n";
    std::cout.flush();
    if ( resid < target_residual )
        return resid;

    for ( size_t k = 0; k < max_iter; ++k )
    {
        Ap = A*p;
        real rz    = dot(r,z);
        real alpha = rz/dot(p,Ap);
        x = x + alpha*p;
        r = r - alpha*Ap;

        norm_r = norm(r);
        resid  = relative ? norm_r/norm_b : norm_r;
        if ( (k + 1) % 10 == 0 )
        {
            std::cout << "CG: Iteration: " << std::setw(4) << k + 1 << ", "
                      << "Residual:  " << std::setw(12) << std::scientific << resid << ".\n";
            std::cout.flush();
        }
        if ( resid <= target_residual )
        {
            std::cout << "CG: Target residual "  << target_residual
                      << " reached after " << k + 1 << " iterations.\n";
            std::cout.flush();
            return resid;
        }

        z = Q*r;
        real beta = dot(r,z)/rz;
        p  = z + beta*p;
    }

    return resid;
}

}

