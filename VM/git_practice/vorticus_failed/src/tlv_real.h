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
#ifndef TLV_REAL_H
#define TLV_REAL_H

#include <array>

#include "types.h"

template <size_t n>
class tlv_real
{
public:
    constexpr  tlv_real() noexcept = default;
    constexpr  tlv_real( const tlv_real&  rhs ) noexcept = default;
    constexpr  tlv_real(       tlv_real&& rhs ) noexcept = default;
    constexpr  tlv_real& operator=( const tlv_real&  rhs ) noexcept = default;
    constexpr  tlv_real& operator=(       tlv_real&& rhs ) noexcept = default;
              ~tlv_real() = default;

    constexpr tlv_real( const real rhs ) noexcept;
    constexpr tlv_real( const real rhs, const std::array<real,n> &drhs ) noexcept;

    void operator+=( const tlv_real &rhs ) noexcept;   
    void operator-=( const tlv_real &rhs ) noexcept;   
    void operator*=( const tlv_real &rhs ) noexcept;   
    void operator/=( const tlv_real &rhs ) noexcept;   

    void operator+=( const real rhs ) noexcept;   
    void operator-=( const real rhs ) noexcept;   
    void operator*=( const real rhs ) noexcept;   
    void operator/=( const real rhs ) noexcept;   

    tlv_real operator+( const tlv_real &rhs ) const noexcept;   
    tlv_real operator-( const tlv_real &rhs ) const noexcept;   
    tlv_real operator*( const tlv_real &rhs ) const noexcept;   
    tlv_real operator/( const tlv_real &rhs ) const noexcept;   

    tlv_real operator+( const real rhs ) const noexcept;   
    tlv_real operator-( const real rhs ) const noexcept;   
    tlv_real operator*( const real rhs ) const noexcept;   
    tlv_real operator/( const real rhs ) const noexcept;   

    tlv_real operator-() const noexcept;

    tlv_real& operator=( const real rhs ) noexcept;

    real               val          {};
    std::array<real,n> derivatives  {};
};

using tl_real = tlv_real<3>;

template <size_t n>
tlv_real<n> operator+( const real lhs, const tlv_real<n> &rhs ) noexcept;

template <size_t n>
tlv_real<n> operator-( const real lhs, const tlv_real<n> &rhs ) noexcept;

template <size_t n>
tlv_real<n> operator*( const real lhs, const tlv_real<n> &rhs ) noexcept;

template <size_t n>
tlv_real<n> operator/( const real lhs, const tlv_real<n> &rhs ) noexcept;


template <size_t n>
tlv_real<n> exp( const tlv_real<n> &arg ) noexcept;

template <size_t n>
tlv_real<n> log( const tlv_real<n> &arg ) noexcept;

template <size_t n>
tlv_real<n> sin( const tlv_real<n> &arg ) noexcept;

template <size_t n>
tlv_real<n> cos( const tlv_real<n> &arg ) noexcept;

template <size_t n>
tlv_real<n> tan( const tlv_real<n> &arg ) noexcept;

template <size_t n>
tlv_real<n> sinh( const tlv_real<n> &arg ) noexcept;

template <size_t n>
tlv_real<n> cosh( const tlv_real<n> &arg ) noexcept;

template <size_t n>
tlv_real<n> tanh( const tlv_real<n> &arg ) noexcept;

template <size_t n>
tlv_real<n> asin( const tlv_real<n> &arg ) noexcept;

template <size_t n>
tlv_real<n> acos( const tlv_real<n> &arg ) noexcept;

template <size_t n>
tlv_real<n> atan( const tlv_real<n> &arg ) noexcept;

template <size_t n>
tlv_real<n> atan2( const tlv_real<n> &y, const tlv_real<n> &x ) noexcept;

template <size_t n>
tlv_real<n> asinh( const tlv_real<n> &arg ) noexcept;

template <size_t n>
tlv_real<n> acosh( const tlv_real<n> &arg ) noexcept;

template <size_t n>
tlv_real<n> atanh( const tlv_real<n> &arg ) noexcept;

template <size_t n>
tlv_real<n> sqrt( const tlv_real<n> &arg ) noexcept;

template <size_t n>
tlv_real<n> abs( const tlv_real<n> &arg ) noexcept;

#include "tlv_real.tpp"
#endif

