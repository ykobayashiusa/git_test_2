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
#ifndef GEOMETRY_TENSOR_H
#define GEOMETRY_TENSOR_H

#include "geometry/point.h"

#include <ostream>
#include <iomanip>
#include <algorithm>

namespace geometry
{

class tensor
{
public:
    constexpr tensor() noexcept = default;
    
    constexpr tensor( real r0, real r1, real r2,
                      real r3, real r4, real r5,
                      real r6, real r7, real r8 ) noexcept;

              real& operator()( size_t i, size_t j )       noexcept;
    constexpr real  operator()( size_t i, size_t j ) const noexcept;

    constexpr tensor operator-() const noexcept;

    constexpr       point operator*( const       point &rhs ) const noexcept;

    constexpr tensor operator*( const tensor &rhs ) const noexcept;
    constexpr tensor operator+( const tensor &rhs ) const noexcept;
    constexpr tensor operator-( const tensor &rhs ) const noexcept;

    constexpr tensor operator/( const real rhs ) const noexcept;
    constexpr tensor operator*( const real rhs ) const noexcept;

    constexpr real   det() const noexcept;
    constexpr tensor adj() const noexcept;
    constexpr tensor inv() const noexcept;

    void operator+=( const tensor &rhs ) noexcept;
    void operator-=( const tensor &rhs ) noexcept;

    void operator/=( const real rhs ) noexcept;
    void operator*=( const real rhs ) noexcept;

    void transpose() noexcept;
    real norm() const;
    
    real data[3][3] { {0,0,0}, {0,0,0}, {0,0,0} };
};

constexpr tensor tensor_colwise( const point &col1, const point &col2, const point &col3 );
constexpr tensor tensor_rowwise( const point &row1, const point &row2, const point &row3 );

constexpr tensor operator*( real lhs, const tensor &rhs ) noexcept { return rhs*lhs; }

constexpr tensor  dyad_prod( const point&  lhs, const point&  rhs ) noexcept;
constexpr tensor cross_prod( const point&  lhs, const tensor& rhs ) noexcept;
constexpr tensor cross_prod( const tensor& lhs, const point&  rhs ) noexcept;

constexpr real double_contraction( const tensor &lhs, const tensor &rhs ) noexcept;

std::ostream& operator<<( std::ostream& str, const tensor& rhs );

///////////////////////
// Member functions. //
///////////////////////

constexpr
tensor::tensor( real r0, real r1, real r2,
                real r3, real r4, real r5,
                real r6, real r7, real r8 ) noexcept:
data  { {r0,r1,r2}, {r3,r4,r5}, {r6,r7,r8} }
{}

inline
real& tensor::operator()( size_t i, size_t j ) noexcept
{
    return data[i][j];
}

constexpr
real tensor::operator()( size_t i, size_t j ) const noexcept
{
    return data[i][j];
}

constexpr
tensor tensor::operator-() const noexcept
{
    return tensor { -data[0][0], -data[0][1], -data[0][2],
                    -data[1][0], -data[1][1], -data[1][2],
                    -data[2][0], -data[2][1], -data[2][2] };
}

constexpr
point tensor::operator*( const point &rhs ) const noexcept
{
    return point
    {
        data[0][0]*rhs.x + data[0][1]*rhs.y + data[0][2]*rhs.z,
        data[1][0]*rhs.x + data[1][1]*rhs.y + data[1][2]*rhs.z,
        data[2][0]*rhs.x + data[2][1]*rhs.y + data[2][2]*rhs.z
    };
}

constexpr
tensor tensor::operator*( const tensor &rhs ) const noexcept
{
    return tensor
    {
            data[0][0]*rhs.data[0][0] + data[0][1]*rhs.data[1][0] + data[0][2]*rhs.data[2][0],
            data[0][0]*rhs.data[0][1] + data[0][1]*rhs.data[1][1] + data[0][2]*rhs.data[2][1],
            data[0][0]*rhs.data[0][2] + data[0][1]*rhs.data[1][2] + data[0][2]*rhs.data[2][2],
        
            data[1][0]*rhs.data[0][0] + data[1][1]*rhs.data[1][0] + data[1][2]*rhs.data[2][0],
            data[1][0]*rhs.data[0][1] + data[1][1]*rhs.data[1][1] + data[1][2]*rhs.data[2][1],
            data[1][0]*rhs.data[0][2] + data[1][1]*rhs.data[1][2] + data[1][2]*rhs.data[2][2],
        
            data[2][0]*rhs.data[0][0] + data[2][1]*rhs.data[1][0] + data[2][2]*rhs.data[2][0],
            data[2][0]*rhs.data[0][1] + data[2][1]*rhs.data[1][1] + data[2][2]*rhs.data[2][1],
            data[2][0]*rhs.data[0][2] + data[2][1]*rhs.data[1][2] + data[2][2]*rhs.data[2][2]
    };
}

inline
void tensor::operator+=( const tensor &rhs ) noexcept
{
    data[0][0] += rhs.data[0][0];
    data[0][1] += rhs.data[0][1];
    data[0][2] += rhs.data[0][2];

    data[1][0] += rhs.data[1][0];
    data[1][1] += rhs.data[1][1];
    data[1][2] += rhs.data[1][2];

    data[2][0] += rhs.data[2][0];
    data[2][1] += rhs.data[2][1];
    data[2][2] += rhs.data[2][2];
}

constexpr
tensor tensor::operator+( const tensor& rhs ) const noexcept
{
    return tensor
    {
            data[0][0] + rhs.data[0][0],
            data[0][1] + rhs.data[0][1],
            data[0][2] + rhs.data[0][2],

            data[1][0] + rhs.data[1][0],
            data[1][1] + rhs.data[1][1],
            data[1][2] + rhs.data[1][2],

            data[2][0] + rhs.data[2][0],
            data[2][1] + rhs.data[2][1],
            data[2][2] + rhs.data[2][2]
    };
}

inline
void tensor::operator-=( const tensor &rhs ) noexcept
{
    data[0][0] -= rhs.data[0][0];
    data[0][1] -= rhs.data[0][1];
    data[0][2] -= rhs.data[0][2];

    data[1][0] -= rhs.data[1][0];
    data[1][1] -= rhs.data[1][1];
    data[1][2] -= rhs.data[1][2];

    data[2][0] -= rhs.data[2][0];
    data[2][1] -= rhs.data[2][1];
    data[2][2] -= rhs.data[2][2];
}

constexpr
tensor tensor::operator-( const tensor& rhs ) const noexcept
{
    return tensor
    {
            data[0][0] - rhs.data[0][0],
            data[0][1] - rhs.data[0][1],
            data[0][2] - rhs.data[0][2],
        
            data[1][0] - rhs.data[1][0],
            data[1][1] - rhs.data[1][1],
            data[1][2] - rhs.data[1][2],

            data[2][0] - rhs.data[2][0],
            data[2][1] - rhs.data[2][1],
            data[2][2] - rhs.data[2][2]
    };
}

constexpr
real tensor::det() const noexcept
{
    return data[0][0]*data[1][1]*data[2][2] +
           data[0][1]*data[1][2]*data[2][0] +
           data[0][2]*data[1][0]*data[2][1] -
           data[2][0]*data[1][1]*data[0][2] -
           data[2][1]*data[1][2]*data[0][0] -
           data[2][2]*data[1][0]*data[0][1];
           
}

constexpr
tensor tensor::adj() const noexcept
{
    return tensor
    {
        data[1][1]*data[2][2] - data[1][2]*data[2][1],
        data[0][2]*data[2][1] - data[0][1]*data[2][2],
        data[0][1]*data[1][2] - data[0][2]*data[1][1],

        data[1][2]*data[2][0] - data[1][0]*data[2][2],
        data[0][0]*data[2][2] - data[0][2]*data[2][0],
        data[0][2]*data[1][0] - data[0][0]*data[1][2],

        data[1][0]*data[2][1] - data[1][1]*data[2][0],
        data[0][1]*data[2][0] - data[0][0]*data[2][1],
        data[0][0]*data[1][1] - data[0][1]*data[1][0]
    };
}

constexpr
tensor tensor::inv() const noexcept
{
    return adj() / det();
}

inline
void tensor::transpose() noexcept
{
    std::swap( data[1][0], data[0][1] );
    std::swap( data[2][0], data[0][2] );
    std::swap( data[2][1], data[1][2] );
}

inline
void tensor::operator/=( const real rhs ) noexcept
{
    data[0][0] /= rhs;
    data[0][1] /= rhs;
    data[0][2] /= rhs;

    data[1][0] /= rhs;
    data[1][1] /= rhs;
    data[1][2] /= rhs;

    data[2][0] /= rhs;
    data[2][1] /= rhs;
    data[2][2] /= rhs;
}

inline
void tensor::operator*=( const real rhs ) noexcept
{
    data[0][0] *= rhs;
    data[0][1] *= rhs;
    data[0][2] *= rhs;

    data[1][0] *= rhs;
    data[1][1] *= rhs;
    data[1][2] *= rhs;

    data[2][0] *= rhs;
    data[2][1] *= rhs;
    data[2][2] *= rhs;
}

constexpr
tensor tensor::operator/( const real rhs ) const noexcept
{
    return tensor
    {
            data[0][0] / rhs,
            data[0][1] / rhs,
            data[0][2] / rhs,

            data[1][0] / rhs,
            data[1][1] / rhs,
            data[1][2] / rhs,

            data[2][0] / rhs,
            data[2][1] / rhs,
            data[2][2] / rhs
    };
}

constexpr
tensor tensor::operator*( const real rhs ) const noexcept
{
    return tensor
    {
            data[0][0] * rhs,
            data[0][1] * rhs,
            data[0][2] * rhs,

            data[1][0] * rhs,
            data[1][1] * rhs,
            data[1][2] * rhs,

            data[2][0] * rhs,
            data[2][1] * rhs,
            data[2][2] * rhs
    };
}

constexpr
tensor tensor_colwise( const point &col1, const point &col2, const point &col3 )
{
    return tensor { col1.x, col2.x, col3.x,
                    col1.y, col2.y, col3.y,
                    col1.z, col2.z, col3.z };
}

constexpr
tensor tensor_rowwise( const point &row1, const point &row2, const point &row3 )
{
    return tensor { row1.x, row1.y, row1.z,
                    row2.x, row2.y, row2.z,
                    row3.x, row3.y, row3.z };
}

constexpr
tensor dyad_prod( const point& a, const point& b ) noexcept
{
    return tensor
    {
        a.x*b.x, a.x*b.y, a.x*b.z,
        a.y*b.x, a.y*b.y, a.y*b.z,
        a.z*b.x, a.z*b.y, a.z*b.z
    };
}

constexpr
tensor cross_prod( const tensor& l, const point& r ) noexcept
{
    return tensor
    (
        l(1,0)*r.z - l(2,0)*r.y,    l(1,1)*r.z - l(2,1)*r.y,    l(1,2)*r.z - l(2,2)*r.y, 
        l(2,0)*r.x - l(0,0)*r.z,    l(2,1)*r.x - l(0,1)*r.z,    l(2,2)*r.x - l(0,2)*r.z,
        l(0,0)*r.y - l(1,0)*r.x,    l(0,1)*r.y - l(1,1)*r.x,    l(0,2)*r.y - l(1,2)*r.x
    );
}

constexpr
tensor cross_prod( const point &l, const tensor &r ) noexcept
{
    return tensor
    (
        r(2,0)*l.y - r(1,0)*l.z,    r(2,1)*l.y - r(1,1)*l.z,    r(2,2)*l.y - r(1,2)*l.z,
        r(0,0)*l.z - r(2,0)*l.x,    r(0,1)*l.z - r(2,1)*l.x,    r(0,2)*l.z - r(2,2)*l.x,
        r(1,0)*l.x - r(0,0)*l.y,    r(1,1)*l.x - r(0,1)*l.y,    r(1,2)*l.x - r(0,2)*l.y
    );
}

constexpr real double_contraction( const tensor &lhs, const tensor &rhs ) noexcept
{
    return lhs.data[0][0]*rhs.data[0][0] + lhs.data[0][1]*rhs.data[0][1] + lhs.data[0][2]*rhs.data[0][2] +
           lhs.data[1][0]*rhs.data[1][0] + lhs.data[1][1]*rhs.data[1][1] + lhs.data[1][2]*rhs.data[1][2] +
           lhs.data[2][0]*rhs.data[2][0] + lhs.data[2][1]*rhs.data[2][1] + lhs.data[2][2]*rhs.data[2][2];
}

inline
std::ostream& operator<<( std::ostream& str, const tensor& rhs )
{
    using std::setw;
    str << setw(20) << rhs.data[0][0] << setw(20) << rhs.data[0][1] << setw(20) << rhs.data[0][2] << '\n'
        << setw(20) << rhs.data[1][0] << setw(20) << rhs.data[1][1] << setw(20) << rhs.data[1][2] << '\n'
        << setw(20) << rhs.data[2][0] << setw(20) << rhs.data[2][1] << setw(20) << rhs.data[2][2] << '\n';
    return str;
}

}

#endif

