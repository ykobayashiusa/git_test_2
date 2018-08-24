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
#ifndef GEOMETRY_ZSORT_H
#define GEOMETRY_ZSORT_H

#include "geometry/point.h"

namespace geometry
{

template <typename iterator>
void zsort( iterator begin, iterator end );

template <typename iterator, typename position_getter>
void zsort( iterator begin, iterator end, position_getter pos );

}

#include "geometry/zsort.tpp"
#endif

