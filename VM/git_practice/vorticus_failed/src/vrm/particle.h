/*
 * Copyright (C) 2015 Matthias Kirchhart
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
#ifndef VRM_PARTICLE_H
#define VRM_PARTICLE_H

#include <array>
#include "geometry/point.h"

namespace vrm
{

template <uint stages>
struct particle
{
    using point = geometry::point;

    point x;
    point G;
 
    std::array<point,stages> k_x;
    std::array<point,stages> k_G;
};

}

#endif

