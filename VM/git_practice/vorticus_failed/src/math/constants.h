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
#ifndef MATH_CONSTANTS_H
#define MATH_CONSTANTS_H

#include <limits>

#include "types.h"

namespace math
{

constexpr real pi         = 3.14159265358979323846;
constexpr real sqrt_two   = 1.41421356237309504880;
constexpr real sqrt_three = 1.73205080756887729352;

constexpr real pifac      = 1./(4.*pi);

constexpr real eps        = std::numeric_limits<real>::epsilon();
}

#endif

