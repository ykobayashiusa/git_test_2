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
#ifndef MISC_STOPWATCH_H
#define MISC_STOPWATCH_H

#include "types.h"
#include <chrono>

namespace misc
{

class stopwatch
{
public:
	void reset();
	real elapsed();

private:
	using clock = std::chrono::high_resolution_clock;
	clock::time_point t0 { clock::now() };
};


inline
void stopwatch::reset()
{
	t0 = clock::now();
}

inline
real stopwatch::elapsed()
{
	using seconds = std::chrono::duration<real,std::ratio<1,1>>;

	auto tnow = clock::now();
	auto duration = std::chrono::duration_cast<seconds>( tnow - t0 );

	return duration.count();
}

}

#endif
