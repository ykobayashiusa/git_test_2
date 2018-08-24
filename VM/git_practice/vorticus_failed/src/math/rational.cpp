/*
 * Copyright (C) 2018 Matthias Kirchhart
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
#include "math/rational.h"

#include <climits>

namespace math
{

void rational::normalise()
{
    // The following loop determines the greatest common divisor of p and q,
    // using Euclid's algorithm. The result is stored in pp. This also works
    // if p == 0 and in that case later normalises q to 1. Using sub(0,p) and
    // sub(0,q) to flip the signs prevents overflows in the case that we have
    // p == INT_MIN or q == INT_MIN.
    int pp { p < 0 ? sub(0,p) : p };
    int qq { q < 0 ? sub(0,q) : q };
    while ( qq != 0 )
    {
        int tmp { pp%qq };
        pp=qq;
        qq=tmp;
    }

    // Simplify the fraction.
    p /= pp;
    q /= pp;

    // Check for negative denominator.
    if ( q < 0 )
    {
        p = sub(0,p); 
        q = sub(0,q);
    }
}

int rational::add( int a, int b )
{
    int c;
    #ifdef __GNUC__ > 5
    if ( __builtin_sadd_overflow(a,b,&c) )
        throw std::range_error { "Integer overflow in addition." };
    #else
    if ( (b > 0) && (a > (INT_MAX - b)) ||
         (b < 0) && (a < (INT_MIN - b)) )
        throw std::range_error { "Integer overflow in addition." };
    c = a + b;
    #endif
    return c;
}

int rational::sub( int a, int b )
{
    int c;
    #ifdef __GNUC__ > 5
    if ( __builtin_ssub_overflow(a,b,&c) )
        throw std::range_error { "Integer overflow in substraction." };
    #else
    if ( (b > 0) && (a < (INT_MIN + b)) ||
         (b < 0) && (a > (INT_MAX + b)) )
        throw std::range_error { "Integer overflow in substraction." };
    c = a - b;
    #endif
    return c;
}

int rational::mul( int a, int b )
{
    int c;
    #ifdef __GNUC__ > 5
    if ( __builtin_smul_overflow(a,b,&c) )
        throw std::range_error { "Integer overflow in multiplication." };
    #else
    if ( ( (a> 0) && (b> 0) && (a>(INT_MAX/b)) ) ||
         ( (a> 0) && (b<=0) && (b<(INT_MIN/a)) ) ||
         ( (a<=0) && (b> 0) && (a<(INT_MIN/b)) ) ||
         ( (a< 0) && (b<=0) && (b<(INT_MAX/a)) ) )
        throw std::range_error { "Integer overflow in multiplication." };
    c = a*b;
    #endif
    return c;
}

}

