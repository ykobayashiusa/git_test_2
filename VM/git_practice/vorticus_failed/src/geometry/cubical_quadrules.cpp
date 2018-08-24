/*
 * Copyright (C) 2017 Matthias Kirchhart
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
#include "geometry/quadrules.h"

namespace
{

const std::pair<real,real> legendre_rule_01[1]
{
    { 0., 2 }
};

const std::pair<real,real> legendre_rule_02[2]
{
    { -0.5773502691896257, 1 },
    {  0.5773502691896257, 1 }
};

const std::pair<real,real> legendre_rule_03[3]
{
    {  0.0000000000000000, 0.8888888888888888 },
    { -0.7745966692414834, 0.5555555555555556 },
    {  0.7745966692414834, 0.5555555555555556 }
};

const std::pair<real,real> legendre_rule_04[4]
{
    { -0.3399810435848563, 0.6521451548625461 },
    {  0.3399810435848563, 0.6521451548625461 },
    { -0.8611363115940526, 0.3478548451374538 },
    {  0.8611363115940526, 0.3478548451374538 }
};

const std::pair<real,real> legendre_rule_05[5]
{
    {  0.0000000000000000, 0.5688888888888889 },
    { -0.5384693101056831, 0.4786286704993665 },
    {  0.5384693101056831, 0.4786286704993665 },
    { -0.9061798459386640, 0.2369268850561891 },
    {  0.9061798459386640, 0.2369268850561891 }
};

const std::pair<real,real> legendre_rule_06[6]
{
    { -0.6612093864662645, 0.3607615730481386 },
    {  0.6612093864662645, 0.3607615730481386 },
    { -0.2386191860831969, 0.4679139345726910 },
    {  0.2386191860831969, 0.4679139345726910 },
    { -0.9324695142031521, 0.1713244923791704 },
    {  0.9324695142031521, 0.1713244923791704 }
};

const std::pair<real,real> legendre_rule_07[7]
{
    {  0.0000000000000000, 0.4179591836734694 },
    { -0.4058451513773972, 0.3818300505051189 },
    {  0.4058451513773972, 0.3818300505051189 },
    { -0.7415311855993945, 0.2797053914892766 },
    {  0.7415311855993945, 0.2797053914892766 },
    { -0.9491079123427585, 0.1294849661688697 },
    {  0.9491079123427585, 0.1294849661688697 }
};

const std::pair<real,real> legendre_rule_08[8]
{
    { -0.1834346424956498, 0.3626837833783620 },
    {  0.1834346424956498, 0.3626837833783620 },
    { -0.5255324099163290, 0.3137066458778873 },
    {  0.5255324099163290, 0.3137066458778873 },
    { -0.7966664774136267, 0.2223810344533745 },
    {  0.7966664774136267, 0.2223810344533745 },
    { -0.9602898564975363, 0.1012285362903763 },
    {  0.9602898564975363, 0.1012285362903763 }
};

const std::pair<real,real> legendre_rule_09[9]
{
    {  0.0000000000000000, 0.3302393550012598 },
    { -0.8360311073266358, 0.1806481606948574 },
    {  0.8360311073266358, 0.1806481606948574 },
    { -0.9681602395076261, 0.0812743883615744 },
    {  0.9681602395076261, 0.0812743883615744 },
    { -0.3242534234038089, 0.3123470770400029 },
    {  0.3242534234038089, 0.3123470770400029 },
    { -0.6133714327005904, 0.2606106964029354 },
    {  0.6133714327005904, 0.2606106964029354 }
};

const std::pair<real,real> legendre_rule_10[10]
{
    { -0.1488743389816312, 0.2955242247147529 },
    {  0.1488743389816312, 0.2955242247147529 },
    { -0.4333953941292472, 0.2692667193099963 },
    {  0.4333953941292472, 0.2692667193099963 },
    { -0.6794095682990244, 0.2190863625159820 },
    {  0.6794095682990244, 0.2190863625159820 },
    { -0.8650633666889845, 0.1494513491505806 },
    {  0.8650633666889845, 0.1494513491505806 },
    { -0.9739065285171717, 0.0666713443086881 },
    {  0.9739065285171717, 0.0666713443086881 }
};

}

namespace geometry
{

const cubical_quadrule get_cubical_quadrule( size_t degree )
{
    size_t N;
    const std::pair<real,real> *rule;

    switch ( degree )
    {
    case  0: rule = legendre_rule_01; N =  1; break;
    case  1: rule = legendre_rule_01; N =  1; break;
    case  2: rule = legendre_rule_02; N =  2; break;
    case  3: rule = legendre_rule_02; N =  2; break;
    case  4: rule = legendre_rule_03; N =  3; break;
    case  5: rule = legendre_rule_03; N =  3; break;
    case  6: rule = legendre_rule_04; N =  4; break;
    case  7: rule = legendre_rule_04; N =  4; break;
    case  8: rule = legendre_rule_05; N =  5; break;
    case  9: rule = legendre_rule_05; N =  5; break;
    case 10: rule = legendre_rule_06; N =  6; break;
    case 11: rule = legendre_rule_06; N =  6; break;
    case 12: rule = legendre_rule_07; N =  7; break;
    case 13: rule = legendre_rule_07; N =  7; break;
    case 14: rule = legendre_rule_08; N =  8; break;
    case 15: rule = legendre_rule_08; N =  8; break;
    case 16: rule = legendre_rule_09; N =  9; break;
    case 17: rule = legendre_rule_09; N =  9; break;
    case 18: rule = legendre_rule_10; N = 10; break;
    case 19: rule = legendre_rule_10; N = 10; break;
    default:
        throw std::logic_error { "Requested too high order cubical quadrature rule." }; 
    }

    cubical_quadrule result( N*N*N );
    for ( size_t i = 0; i < N; ++i )
    for ( size_t j = 0; j < N; ++j )
    for ( size_t k = 0; k < N; ++k )
    {
        result[ i*N*N + j*N + k ]= { point { rule[i].first,   rule[j].first,   rule[k].first  },
                                     real  { rule[i].second * rule[j].second * rule[k].second } };
    }
    return std::move(result);
}

cubical_quadrule get_refined_cubical_quadrule( size_t degree, size_t level )
{
    size_t N;
    const std::pair<real,real> *rule;

    switch ( degree )
    {
    case  0: rule = legendre_rule_01; N =  1; break;
    case  1: rule = legendre_rule_01; N =  1; break;
    case  2: rule = legendre_rule_02; N =  2; break;
    case  3: rule = legendre_rule_02; N =  2; break;
    case  4: rule = legendre_rule_03; N =  3; break;
    case  5: rule = legendre_rule_03; N =  3; break;
    case  6: rule = legendre_rule_04; N =  4; break;
    case  7: rule = legendre_rule_04; N =  4; break;
    case  8: rule = legendre_rule_05; N =  5; break;
    case  9: rule = legendre_rule_05; N =  5; break;
    case 10: rule = legendre_rule_06; N =  6; break;
    case 11: rule = legendre_rule_06; N =  6; break;
    case 12: rule = legendre_rule_07; N =  7; break;
    case 13: rule = legendre_rule_07; N =  7; break;
    case 14: rule = legendre_rule_08; N =  8; break;
    case 15: rule = legendre_rule_08; N =  8; break;
    case 16: rule = legendre_rule_09; N =  9; break;
    case 17: rule = legendre_rule_09; N =  9; break;
    case 18: rule = legendre_rule_10; N = 10; break;
    case 19: rule = legendre_rule_10; N = 10; break;
    default:
        throw std::logic_error { "Requested too high order cubical quadrature rule." }; 
    }

    size_t M   { 1ul << level };
    real   dx  { 2./M };
    std::vector< std::pair<real,real> > rule_1d( N*M );
    for ( size_t i = 0; i < M; ++i )
    {
        real a = -1 + (i  )*dx;
        real b = -1 + (i+1)*dx;
        for ( size_t j = 0; j < N; ++j )
        {
            rule_1d[ N*i + j ].first  = (dx/2)*rule[ j ].first + (a+b)/2;
            rule_1d[ N*i + j ].second = (dx/2)*rule[ j ].second;
        }
    }

    N = N*M;
    cubical_quadrule result( N*N*N );
    for ( size_t i = 0; i < N; ++i )
    for ( size_t j = 0; j < N; ++j )
    for ( size_t k = 0; k < N; ++k )
    {
        result[ i*N*N + j*N + k ]= { point { rule_1d[i].first,   rule_1d[j].first,   rule_1d[k].first  },
                                     real  { rule_1d[i].second * rule_1d[j].second * rule_1d[k].second } };
    }
    return std::move(result);
}

}

