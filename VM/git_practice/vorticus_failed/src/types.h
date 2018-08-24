/*
 * Copyright (C) 2013 Matthias Kirchhart
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

#ifndef TYPES_H
#define TYPES_H

#include <cmath>
#include <limits>
#include <complex>
#include <stdexcept>


using uint  = unsigned int;
using uchar = unsigned char;

using real  = double; // Should be defined as a 64-bit float.
using cmplx = std::complex<real>;

#include "tlv_real.h"

template <typename target, typename source>
target narrow_cast( source val )
{
    target res { static_cast<target>(val) };
    if ( static_cast<source>(res) != val )
        throw std::runtime_error { "narrow_cast<>() failed!" };
    return res;
}

template <typename target, typename source>
bool narrow_castable( source val )
{
    target res { static_cast<target>(val) };
    return static_cast<source>(res) == val;
}

#endif

