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
#ifndef MATH_RATIONAL_H
#define MATH_RATIONAL_H

#include <stdexcept>

namespace math
{

class rational
{
public:
    constexpr rational() noexcept = default;
    constexpr rational( int pp ) noexcept;
              rational( int pp, int qq );
    constexpr rational( const rational &&rhs ) noexcept = default;
    constexpr rational(       rational  &rhs ) noexcept = default;

    rational& operator=( const rational &&rhs ) noexcept = default;
    rational& operator=( const rational  &rhs ) noexcept = default;
    rational& operator=( int rhs ) noexcept;

    constexpr bool operator==( const rational &rhs ) const noexcept;
    constexpr bool operator!=( const rational &rhs ) const noexcept;

    bool operator< ( const rational &rhs ) const;
    bool operator> ( const rational &rhs ) const;
    bool operator<=( const rational &rhs ) const;
    bool operator>=( const rational &rhs ) const;

    rational& operator+=( const rational &rhs );
    rational& operator-=( const rational &rhs );
    rational& operator*=( const rational &rhs );
    rational& operator/=( const rational &rhs );

    rational& operator+=( int rhs );
    rational& operator-=( int rhs );
    rational& operator*=( int rhs );
    rational& operator/=( int rhs );

    rational operator+( const rational &rhs ) const;
    rational operator-( const rational &rhs ) const;
    rational operator*( const rational &rhs ) const;
    rational operator/( const rational &rhs ) const;

    rational operator+( int rhs ) const;
    rational operator-( int rhs ) const;
    rational operator*( int rhs ) const;
    rational operator/( int rhs ) const;

    rational operator-() const;
    rational operator+() const noexcept;

private:
    void normalise() noexcept;

    static int add(int,int);
    static int sub(int,int);
    static int mul(int,int);

private:
    int p {0};
    int q {1};
};

}

#include "math/rational.tpp"
#endif

