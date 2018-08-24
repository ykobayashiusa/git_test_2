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
#include <limits>
#include <iterator>

namespace geometry
{

template <typename iterator>
void zsort( iterator begin, iterator end )
{
    if ( end - begin <= 1 ) return;

    constexpr real MAX { std::numeric_limits<real>::max() };
    point min_corner {  MAX,  MAX,  MAX };
    point max_corner { -MAX, -MAX, -MAX };

    for ( iterator i = begin; i != end; ++i )
    {
        min_corner = min_components( *i, min_corner );
        max_corner = max_components( *i, max_corner ); 
    }

    const point centre = (min_corner + max_corner)/2;
    
    using reference = typename std::iterator_traits<iterator>::reference;
    auto xless = [centre]( reference i ) { return i.x < centre.x; };
    auto yless = [centre]( reference i ) { return i.y < centre.y; };
    auto zless = [centre]( reference i ) { return i.z < centre.z; };
    
    iterator oct[9];
    oct[0] = begin;
    oct[8] = end;
        
    oct[4] = std::partition( oct[0], oct[8], xless );
    oct[2] = std::partition( oct[0], oct[4], yless );
    oct[6] = std::partition( oct[4], oct[8], yless );
    oct[1] = std::partition( oct[0], oct[2], zless );
    oct[3] = std::partition( oct[2], oct[4], zless );
    oct[5] = std::partition( oct[4], oct[6], zless );
    oct[7] = std::partition( oct[6], oct[8], zless );

    #pragma omp parallel for
    for ( size_t i = 0; i < 8; ++i )
        zsort( oct[i], oct[i+1] );
}

template <typename iterator, typename position_getter>
void zsort( iterator begin, iterator end, position_getter pos )
{
    if ( end - begin <= 1 ) return;

    constexpr real MAX { std::numeric_limits<real>::max() };
    point min_corner {  MAX,  MAX,  MAX };
    point max_corner { -MAX, -MAX, -MAX };

    for ( iterator i = begin; i != end; ++i )
    {
        min_corner = min_components( pos(*i), min_corner );
        max_corner = max_components( pos(*i), max_corner ); 
    }

    const point centre = (min_corner + max_corner)/2;
    
    using reference = typename std::iterator_traits<iterator>::reference;
    auto xless = [centre,pos]( reference i ) { return pos(i).x < centre.x; };
    auto yless = [centre,pos]( reference i ) { return pos(i).y < centre.y; };
    auto zless = [centre,pos]( reference i ) { return pos(i).z < centre.z; };
    
    iterator oct[9];
    oct[0] = begin;
    oct[8] = end;
        
    oct[4] = std::partition( oct[0], oct[8], xless );
    oct[2] = std::partition( oct[0], oct[4], yless );
    oct[6] = std::partition( oct[4], oct[8], yless );
    oct[1] = std::partition( oct[0], oct[2], zless );
    oct[3] = std::partition( oct[2], oct[4], zless );
    oct[5] = std::partition( oct[4], oct[6], zless );
    oct[7] = std::partition( oct[6], oct[8], zless );

    #pragma omp parallel for
    for ( size_t i = 0; i < 8; ++i )
        zsort( oct[i], oct[i+1], pos );
}

}

