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
#ifndef MISC_RANDOM_H
#define MISC_RANDOM_H

#include "types.h"

#include <random>
#include <functional>

namespace misc
{

class random_real
{
public:
	random_real( real min, real max );

	real operator()() const;
private:
	std::function<real()> r;
};




inline
random_real::random_real( real min, real max ):
 r( std::bind( std::uniform_real_distribution<>(min,max),
               std::default_random_engine() ) )
{}

inline
real random_real::operator()() const
{
	return r();
}

}

#endif

