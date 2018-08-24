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

namespace math
{

constexpr
rational::rational( int pp ) noexcept:
p { pp }, q { 1 }
{}

rational::rational( int pp, int qq ):
p { pp }, q { qq }
{
    if ( q == 0 )
        throw std::domain_error { "Division by zero." };

    normalise();
}

inline
rational& rational::operator=( int rhs ) noexcept
{
    p = rhs;
    q = 1;
}

bool rational::operator==( const rational &rhs ) const noexcept
{
    return (p==rhs.p) && (q==rhs.q);
}

constexpr
bool rational::operator!=( const rational &rhs ) const noexcept
{
    return (p!=rhs.p) || (q!=rhs.q);
}

inline
rational& rational::operator+=( const rational &rhs )
{
    p = add(mul(p,rhs.q),mul(rhs.p,q));
    q = mul(q,rhs.q);
    normalise();
    return *this;
}

inline
rational& rational::operator-=( const rational &rhs )
{
    p = sub(mul(p,rhs.q),mul(rhs.p,q));
    q = mul(q,rhs.q);
    normalise();
    return *this;
}

inline
rational& rational::operator*=( const rational &rhs )
{
    p = mul(p,rhs.p);
    q = mul(q,rhs.q);
    normalise();
    return *this;
}

inline
rational& rational::operator/=( const rational &rhs )
{
    if ( rhs.p == 0 )
        throw std::domain_error { "Division by zero." };

    if ( this != &rhs )
    {
        p = mul(p,rhs.q);
        q = mul(q,rhs.p);
        normalise();
    }
    else
    {
        p = q = 1;
    }
    return *this;
}

}

