/*
 * Copyright (C) 2017 Matthias Kirchhart
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
#include <iostream>
#include "fem/potentials.h"

#include "math/constants.h"
#include "geometry/quadrules.h"
#include "geometry/tensor.h"

using namespace fem;

int main()
{
    constexpr geometry::tensor id { 1, 0, 0,
                                    0, 1, 0,
                                    0, 0, 1 };

    std::array<point,3> X {{ point {   2,0,0},
                             point {-1.5,2,0},
                             point {  -1,0,0} }};
   
       point   P   { X[2] + point { 0, 0, 1 } };
    tl_point tlP   { P, id };

    real   const_singlet      {}; point grad_const_singlet {};
    real   const_doublet      {}; point grad_const_doublet {};
    real  linear_singlet[3] = {}; point grad_linear_singlet[3] = {};
    real  linear_doublet[3] = {}; point grad_linear_doublet[3] = {};

    point N = cross_prod( X[1]-X[0], X[2]-X[1] );
    real  A = N.r()/2; N /= 2*A;
    triangular_quadrule rule { get_triangular_quadrule(20) };
    for ( size_t i = 0; i < rule.size(); ++i )
    {
        bary2d b { rule[i].b   };
        real   w { rule[i].w*A*math::pifac };
        
        point y = b.z0*X[0]+b.z1*X[1]+b.z2*X[2];
        point R = y-P; real r = R.r(); real r3 = r*r*r; real r5 = r3*r*r;

         const_singlet    += w        / r;
        linear_singlet[0] += w * b.z0 / r; 
        linear_singlet[1] += w * b.z1 / r;
        linear_singlet[2] += w * b.z2 / r;

         const_doublet    += w *        scal_prod(-R,N)/ r3;
        linear_doublet[0] += w * b.z0 * scal_prod(-R,N) / r3;
        linear_doublet[1] += w * b.z1 * scal_prod(-R,N) / r3;
        linear_doublet[2] += w * b.z2 * scal_prod(-R,N) / r3;

         grad_const_singlet    += w *        R / r3;
        grad_linear_singlet[0] += w * b.z0 * R / r3;
        grad_linear_singlet[1] += w * b.z1 * R / r3;
        grad_linear_singlet[2] += w * b.z2 * R / r3;

        geometry::tensor T = ( id/r3 - 3*dyad_prod(R,R)/r5);
          grad_const_doublet    += w *        T*N;
         grad_linear_doublet[0] += w * b.z0 * T*N;
         grad_linear_doublet[1] += w * b.z1 * T*N;
         grad_linear_doublet[2] += w * b.z2 * T*N;
    }

    std::cout << "Potential evaluation: " << std::endl;
    auto const_analytic = layer_potentials<0>( X, P );
    std::cout << "Quadrature result const singlet: " << const_singlet            << std::endl;
    std::cout << "Analytic   result const singlet: " << const_analytic.first[0]  << std::endl;
    std::cout << "Quadrature result const doublet: " << const_doublet            << std::endl;
    std::cout << "Analytic   result const doublet: " << const_analytic.second[0] << std::endl;

    auto linear_analytic = layer_potentials<1>( X, P );
    std::cout << "Quadrature result linear singlet: "
              << linear_singlet[0] << " "
              << linear_singlet[1] << " "
              << linear_singlet[2] << ".\n";
    std::cout << "Analytic   result linear singlet: "
              << linear_analytic.first[0] << " " 
              << linear_analytic.first[1] << " "
              << linear_analytic.first[2] << ".\n";
    
    std::cout << "Quadrature result linear doublet: " 
              << linear_doublet[0]  << " "
              << linear_doublet[1]  << " "
              << linear_doublet[2]  << ".\n";
    std::cout << "Analytic   result linear doublet: "
              << linear_analytic.second[0] << " "
              << linear_analytic.second[1] << " "
              << linear_analytic.second[2] << ".\n";

    std::cout << "=====================================\n";
    std::cout << "Gradient evaluation: \n";

    auto const_grad_analytic = layer_potentials<0>( X, tlP );
    std::cout << "Quadrature result const singlet: " << grad_const_singlet << std::endl;
    std::cout << "Analytic   result const singlet: " << grad(const_grad_analytic.first[0]) << std::endl;
    std::cout << "Quadrature result const doublet: " << grad_const_doublet << std::endl;
    std::cout << "Analytic   result const doublet: " << grad(const_grad_analytic.second[0]) << std::endl;
    
    auto linear_grad_analytic = layer_potentials<1>( X, tlP );
    std::cout << "Quadrature result linear singlet: \n"
              << tensor_rowwise( grad_linear_singlet[0], grad_linear_singlet[1], grad_linear_singlet[2] ) << std::endl;
    std::cout << "Analytic   result linear singlet: \n"
              << tensor_rowwise( grad(linear_grad_analytic.first[0]),
                                 grad(linear_grad_analytic.first[1]),
                                 grad(linear_grad_analytic.first[2]) ) << std::endl;
    
    std::cout << "Quadrature result linear doublet: \n"
              << tensor_rowwise( grad_linear_doublet[0], grad_linear_doublet[1], grad_linear_doublet[2] ) << std::endl;
    std::cout << "Analytic   result linear doublet: \n"
              << tensor_rowwise( grad(linear_grad_analytic.second[0]),
                                 grad(linear_grad_analytic.second[1]),
                                 grad(linear_grad_analytic.second[2]) ) << std::endl;
    
}

