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

template <size_t n> constexpr
tlv_real<n>::tlv_real( const real rhs ) noexcept:
val { rhs }
{}

template <size_t n> constexpr
tlv_real<n>::tlv_real( const real rhs, const std::array<real,n> &drhs ) noexcept:
val { rhs }, derivatives { drhs }
{}

template <size_t n> inline
void tlv_real<n>::operator+=( const tlv_real &rhs ) noexcept
{
    for ( size_t i = 0; i < n; ++i )
    {
        derivatives[ i ] += rhs.derivatives[ i ];
    }
    val += rhs.val;
}

template <size_t n> inline 
void tlv_real<n>::operator-=( const tlv_real &rhs ) noexcept
{
    for ( size_t i = 0; i < n; ++i )
    {
        derivatives[ i ] -= rhs.derivatives[ i ];
    }
    val -= rhs.val;
}

template <size_t n> inline
void tlv_real<n>::operator*=( const tlv_real &rhs ) noexcept
{
    for ( size_t i = 0; i < n; ++i )
    {
        derivatives[ i ] = derivatives[ i ]*rhs.val + val*rhs.derivatives[ i ];
    }
    val *= rhs.val;
}

template <size_t n> inline
void tlv_real<n>::operator/=( const tlv_real &rhs ) noexcept
{
    // To avoid too many expensive divisions, we store the reciprocal of rhs.val.
    real q  = 1/rhs.val;
    real qq = q*q;
    for ( size_t i = 0; i < n; ++i )
    {
        derivatives[ i ] = qq*(derivatives[ i ]*rhs.val - val*rhs.derivatives[ i ]);
    }
    val *= q;
}

template <size_t n> inline
void tlv_real<n>::operator+=( const real rhs ) noexcept
{
    val += rhs;
}

template <size_t n> inline
void tlv_real<n>::operator-=( const real rhs ) noexcept
{
    val -= rhs;
}

template <size_t n> inline
void tlv_real<n>::operator*=( const real rhs ) noexcept
{
    for ( size_t i = 0; i < n; ++i )
    {
        derivatives[ i ] *= rhs;
    }
    val *= rhs;
}

template <size_t n> inline
void tlv_real<n>::operator/=( const real rhs ) noexcept
{
    real q = 1/rhs;
    for ( size_t i = 0; i < n; ++i )
    {
        derivatives[ i ] *= q;
    }
    val *= q;
}

template <size_t n> inline
tlv_real<n> tlv_real<n>::operator+( const tlv_real &rhs ) const noexcept
{
    tlv_real<n> result { *this };
    result += rhs;
    return result;
}

template <size_t n> inline
tlv_real<n> tlv_real<n>::operator-( const tlv_real &rhs ) const noexcept
{
    tlv_real<n> result { *this };
    result -= rhs;
    return result;
}

template <size_t n> inline
tlv_real<n> tlv_real<n>::operator*( const tlv_real &rhs ) const noexcept
{
    tlv_real<n> result { *this };
    result *= rhs;
    return result;
}

template <size_t n> inline
tlv_real<n> tlv_real<n>::operator/( const tlv_real &rhs ) const noexcept
{
    tlv_real<n> result { *this };
    result /= rhs;
    return result;
}

template <size_t n> inline
tlv_real<n> tlv_real<n>::operator+( const real rhs ) const noexcept
{
    tlv_real<n> result { *this };
    result += rhs;
    return result;
}

template <size_t n> inline
tlv_real<n> tlv_real<n>::operator-( const real rhs ) const noexcept
{
    tlv_real<n> result { *this };
    result -= rhs;
    return result;
}

template <size_t n> inline
tlv_real<n> tlv_real<n>::operator*( const real rhs ) const noexcept
{
    tlv_real<n> result { *this };
    result *= rhs;
    return result;
}

template <size_t n> inline
tlv_real<n> tlv_real<n>::operator/( const real rhs ) const noexcept
{
    tlv_real<n> result { *this };
    result /= rhs;
    return result;
}

template <size_t n> inline
tlv_real<n> tlv_real<n>::operator-() const noexcept
{
    tlv_real<n> result;
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] = - derivatives[ i ];
    }
    result.val = - val;
    return result;
}

template <size_t n> inline
tlv_real<n>& tlv_real<n>::operator=( const real rhs ) noexcept
{
    for ( size_t i = 0; i < n; ++i )
    {
        derivatives[ i ] = 0;
    }
    val = rhs;
    return *this;
}

template <size_t n> inline
tlv_real<n> operator+( const real lhs, const tlv_real<n> &rhs ) noexcept
{
    return rhs + lhs;
}

template <size_t n> inline
tlv_real<n> operator-( const real lhs, const tlv_real<n> &rhs ) noexcept
{
    tlv_real<n> result = -rhs;
    result += lhs;
    return result;
}

template <size_t n> inline
tlv_real<n> operator*( const real lhs, const tlv_real<n> &rhs ) noexcept
{
    return rhs*lhs;
}

template <size_t n> inline
tlv_real<n> operator/( const real lhs, const tlv_real<n> &rhs ) noexcept
{
    tlv_real<n> result = lhs;
    result /= rhs;
    return result;
}

template <size_t n> inline
tlv_real<n> exp( const tlv_real<n> &arg ) noexcept
{
    tlv_real<n> result { std::exp(arg.val), arg.derivatives };
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] *= result.val;
    }
    return result;
}

template <size_t n> inline
tlv_real<n> log( const tlv_real<n> &arg ) noexcept
{
    tlv_real<n> result { std::log(arg.val), arg.derivatives };
    real q = 1/arg.val;
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] *= q;
    }
    return result;
}

template <size_t n> inline
tlv_real<n> sin( const tlv_real<n> &arg ) noexcept
{
    tlv_real<n> result { std::sin(arg.val), arg.derivatives };
    const real d  = std::cos(arg.val);
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] *= d;
    }
    return result;
}

template <size_t n> inline
tlv_real<n> cos( const tlv_real<n> &arg ) noexcept
{
    tlv_real<n> result { std::cos(arg.val), arg.derivatives };
    const real d = -std::sin(arg.val);
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] *= d;
    }
    return result;
}

template <size_t n> inline
tlv_real<n> tan( const tlv_real<n> &arg ) noexcept
{
    tlv_real<n> result { std::tan(arg.val), arg.derivatives };
    const real d = 1 + result.val*result.val;
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] *= d;
    }
    return result;
}

template <size_t n> inline
tlv_real<n> sinh( const tlv_real<n> &arg ) noexcept
{
    tlv_real<n> result { std::sinh(arg.val), arg.derivatives };
    //const real d  = std::cosh(arg.val);
    const real d  = std::sqrt( 1 + result.val*result.val );
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] *= d;
    }
    return result;
}

template <size_t n> inline
tlv_real<n> cosh( const tlv_real<n> &arg ) noexcept
{
    tlv_real<n> result { std::cosh(arg.val), arg.derivatives };
    const real d = std::sinh(arg.val);
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] *= d;
    }
    return result;
}

template <size_t n> inline
tlv_real<n> tanh( const tlv_real<n> &arg ) noexcept
{
    tlv_real<n> result { std::tanh(arg.val), arg.derivatives };
    const real d = 1 - result.val*result.val;
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] *= d;
    }
    return result;
}

template <size_t n> inline
tlv_real<n> asin( const tlv_real<n> &arg ) noexcept
{
    tlv_real<n> result { std::asin(arg.val), arg.derivatives };
    const real d  = 1./std::sqrt(1. - arg.val*arg.val);
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] *= d;
    }
    return result;
}

template <size_t n> inline
tlv_real<n> acos( const tlv_real<n> &arg ) noexcept
{
    tlv_real<n> result { std::acos(arg.val), arg.derivatives };
    const real d  = - 1/std::sqrt(1 - arg.val*arg.val);
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] *= d;
    }
    return result;
}

template <size_t n> inline
tlv_real<n> atan( const tlv_real<n> &arg ) noexcept
{
    tlv_real<n> result { std::atan(arg.val), arg.derivatives };
    const real d = 1/(1 + arg.val*arg.val);
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] *= d;
    }
    return result;
}

template <size_t n>
tlv_real<n> atan2( const tlv_real<n> &y, const tlv_real<n> &x ) noexcept
{
    tlv_real<n> result { std::atan2(y.val,x.val) }; 
    real r  = 1./(x.val*x.val + y.val*y.val);
    real dy =  x.val*r; 
    real dx = -y.val*r;
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] = dy*y.derivatives[ i ] +
                                  dx*x.derivatives[ i ];
    }
    return result;
}

template <size_t n> inline
tlv_real<n> asinh( const tlv_real<n> &arg ) noexcept
{
    tlv_real<n> result { std::asinh(arg.val), arg.derivatives };
    const real d = 1/std::sqrt(1 + arg.val*arg.val);
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] *= d;
    }
    return result;
}

template <size_t n> inline
tlv_real<n> acosh( const tlv_real<n> &arg ) noexcept
{
    tlv_real<n> result { std::acosh(arg.val), arg.derivatives };
    const real d = 1/std::sqrt(arg.val*arg.val-1);
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] *= d;
    }
    return result;
}

template <size_t n> inline
tlv_real<n> atanh( const tlv_real<n> &arg ) noexcept
{
    tlv_real<n> result { std::atanh(arg.val), arg.derivatives };
    const real d = 1/(1-arg.val*arg.val);
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] *= d;
    }
    return result;
}

template <size_t n> inline
tlv_real<n> sqrt( const tlv_real<n> &arg ) noexcept
{
    tlv_real<n> result { std::sqrt(arg.val), arg.derivatives };
    const real d = 1/(2*result.val);
    for ( size_t i = 0; i < n; ++i )
    {
        result.derivatives[ i ] *= d;
    }
    return result;
}

template <size_t n>
tlv_real<n> abs( const tlv_real<n> &arg ) noexcept
{
    tlv_real<n> result { std::abs(arg.val), arg.derivatives };
    if ( arg.val < 0 ) 
    {
        for ( size_t i = 0; i < n; ++i )
        {
            result.derivatives[ i ] = - result.derivatives[ i ];
        }
    }
    else if ( arg.val == 0 )
    {
        for ( size_t i = 0; i < n; ++i )
        {
            result.derivatives[ i ] = 0;
        }
    }
    return result;
}

