/*
 * Copyright (C) 2016 Matthias Kirchhart
 *
 * This file is part of vorticus.
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
#ifndef FEM_REFINEMENT_RULES_H
#define FEM_REFINEMENT_RULES_H

#include <array>

namespace fem
{

constexpr std::array<std::array<std::array<unsigned char,2>,4>,94> child_tetrahedra
{{
    {{ {{0,0}}, {{0,1}}, {{0,2}}, {{0,3}} }}, //  0
    {{ {{0,1}}, {{1,1}}, {{1,2}}, {{1,3}} }}, //  1
    {{ {{0,1}}, {{0,2}}, {{1,2}}, {{1,3}} }}, //  2
    {{ {{0,1}}, {{0,2}}, {{0,3}}, {{1,3}} }}, //  3
    {{ {{0,2}}, {{1,2}}, {{2,2}}, {{2,3}} }}, //  4
    {{ {{0,2}}, {{1,2}}, {{1,3}}, {{2,3}} }}, //  5
    {{ {{0,2}}, {{0,3}}, {{1,3}}, {{2,3}} }}, //  6
    {{ {{0,3}}, {{1,3}}, {{2,3}}, {{3,3}} }}, //  7
    {{ {{0,0}}, {{1,1}}, {{2,2}}, {{3,3}} }}, //  8
    {{ {{0,0}}, {{1,1}}, {{2,2}}, {{0,3}} }}, //  9
    {{ {{0,0}}, {{1,1}}, {{2,2}}, {{1,3}} }}, // 10
    {{ {{0,0}}, {{1,1}}, {{2,2}}, {{2,3}} }}, // 11
    {{ {{0,0}}, {{1,1}}, {{0,2}}, {{3,3}} }}, // 12
    {{ {{0,0}}, {{1,1}}, {{0,2}}, {{0,3}} }}, // 13
    {{ {{0,0}}, {{1,1}}, {{0,2}}, {{1,3}} }}, // 14
    {{ {{0,0}}, {{1,1}}, {{1,2}}, {{3,3}} }}, // 15
    {{ {{0,0}}, {{1,1}}, {{1,2}}, {{0,3}} }}, // 16
    {{ {{0,0}}, {{1,1}}, {{1,2}}, {{1,3}} }}, // 17
    {{ {{0,0}}, {{1,1}}, {{2,3}}, {{3,3}} }}, // 18
    {{ {{0,0}}, {{2,2}}, {{1,3}}, {{3,3}} }}, // 19
    {{ {{0,0}}, {{2,2}}, {{1,3}}, {{2,3}} }}, // 20
    {{ {{0,0}}, {{0,1}}, {{2,2}}, {{3,3}} }}, // 21
    {{ {{0,0}}, {{0,1}}, {{2,2}}, {{0,3}} }}, // 22
    {{ {{0,0}}, {{0,1}}, {{2,2}}, {{2,3}} }}, // 23
    {{ {{0,0}}, {{0,1}}, {{0,2}}, {{3,3}} }}, // 24
    {{ {{0,0}}, {{0,1}}, {{2,3}}, {{3,3}} }}, // 25
    {{ {{0,0}}, {{0,2}}, {{1,3}}, {{3,3}} }}, // 26
    {{ {{0,0}}, {{1,2}}, {{2,2}}, {{3,3}} }}, // 27
    {{ {{0,0}}, {{1,2}}, {{2,2}}, {{0,3}} }}, // 28
    {{ {{0,0}}, {{1,2}}, {{2,2}}, {{2,3}} }}, // 29
    {{ {{0,0}}, {{1,2}}, {{1,3}}, {{3,3}} }}, // 30
    {{ {{0,0}}, {{1,2}}, {{1,3}}, {{2,3}} }}, // 31
    {{ {{0,0}}, {{1,2}}, {{2,3}}, {{3,3}} }}, // 32
    {{ {{0,0}}, {{1,3}}, {{2,3}}, {{3,3}} }}, // 33
    {{ {{1,1}}, {{2,2}}, {{0,3}}, {{3,3}} }}, // 34
    {{ {{1,1}}, {{2,2}}, {{0,3}}, {{1,3}} }}, // 35
    {{ {{1,1}}, {{2,2}}, {{0,3}}, {{2,3}} }}, // 36
    {{ {{1,1}}, {{0,2}}, {{2,2}}, {{3,3}} }}, // 37
    {{ {{1,1}}, {{0,2}}, {{2,2}}, {{1,3}} }}, // 38
    {{ {{1,1}}, {{0,2}}, {{2,2}}, {{2,3}} }}, // 39
    {{ {{1,1}}, {{0,2}}, {{1,2}}, {{3,3}} }}, // 40
    {{ {{1,1}}, {{0,2}}, {{1,2}}, {{0,3}} }}, // 41
    {{ {{1,1}}, {{0,2}}, {{1,2}}, {{1,3}} }}, // 42
    {{ {{1,1}}, {{0,2}}, {{0,3}}, {{3,3}} }}, // 43
    {{ {{1,1}}, {{0,2}}, {{0,3}}, {{1,3}} }}, // 44
    {{ {{1,1}}, {{0,2}}, {{0,3}}, {{2,3}} }}, // 45
    {{ {{1,1}}, {{0,2}}, {{2,3}}, {{3,3}} }}, // 46
    {{ {{1,1}}, {{1,2}}, {{0,3}}, {{3,3}} }}, // 47
    {{ {{1,1}}, {{1,2}}, {{0,3}}, {{1,3}} }}, // 48
    {{ {{1,1}}, {{0,3}}, {{2,3}}, {{3,3}} }}, // 49
    {{ {{2,2}}, {{0,3}}, {{1,3}}, {{3,3}} }}, // 50
    {{ {{2,2}}, {{0,3}}, {{1,3}}, {{2,3}} }}, // 51
    {{ {{0,1}}, {{1,1}}, {{2,2}}, {{3,3}} }}, // 52
    {{ {{0,1}}, {{1,1}}, {{2,2}}, {{1,3}} }}, // 53
    {{ {{0,1}}, {{1,1}}, {{2,2}}, {{2,3}} }}, // 54
    {{ {{0,1}}, {{1,1}}, {{1,2}}, {{3,3}} }}, // 55
    {{ {{0,1}}, {{1,1}}, {{2,3}}, {{3,3}} }}, // 56
    {{ {{0,1}}, {{2,2}}, {{0,3}}, {{3,3}} }}, // 57
    {{ {{0,1}}, {{2,2}}, {{0,3}}, {{1,3}} }}, // 58
    {{ {{0,1}}, {{2,2}}, {{0,3}}, {{2,3}} }}, // 59
    {{ {{0,1}}, {{2,2}}, {{1,3}}, {{3,3}} }}, // 60
    {{ {{0,1}}, {{2,2}}, {{1,3}}, {{2,3}} }}, // 61
    {{ {{0,1}}, {{0,2}}, {{2,2}}, {{3,3}} }}, // 62
    {{ {{0,1}}, {{0,2}}, {{2,2}}, {{1,3}} }}, // 63
    {{ {{0,1}}, {{0,2}}, {{2,2}}, {{2,3}} }}, // 64
    {{ {{0,1}}, {{0,2}}, {{1,2}}, {{3,3}} }}, // 65
    {{ {{0,1}}, {{0,2}}, {{1,2}}, {{2,3}} }}, // 66
    {{ {{0,1}}, {{0,2}}, {{0,3}}, {{3,3}} }}, // 67
    {{ {{0,1}}, {{0,2}}, {{0,3}}, {{2,3}} }}, // 68
    {{ {{0,1}}, {{0,2}}, {{1,3}}, {{3,3}} }}, // 69
    {{ {{0,1}}, {{0,2}}, {{2,3}}, {{3,3}} }}, // 70
    {{ {{0,1}}, {{1,2}}, {{2,2}}, {{3,3}} }}, // 71
    {{ {{0,1}}, {{1,2}}, {{2,2}}, {{0,3}} }}, // 72
    {{ {{0,1}}, {{1,2}}, {{2,2}}, {{2,3}} }}, // 73
    {{ {{0,1}}, {{1,2}}, {{0,3}}, {{1,3}} }}, // 74
    {{ {{0,1}}, {{1,2}}, {{1,3}}, {{3,3}} }}, // 75
    {{ {{0,1}}, {{1,2}}, {{1,3}}, {{2,3}} }}, // 76
    {{ {{0,1}}, {{1,2}}, {{2,3}}, {{3,3}} }}, // 77
    {{ {{0,1}}, {{0,3}}, {{1,3}}, {{2,3}} }}, // 78
    {{ {{0,1}}, {{0,3}}, {{2,3}}, {{3,3}} }}, // 79
    {{ {{0,1}}, {{1,3}}, {{2,3}}, {{3,3}} }}, // 80
    {{ {{0,2}}, {{2,2}}, {{1,3}}, {{3,3}} }}, // 81
    {{ {{0,2}}, {{2,2}}, {{1,3}}, {{2,3}} }}, // 82
    {{ {{0,2}}, {{1,2}}, {{2,2}}, {{3,3}} }}, // 83
    {{ {{0,2}}, {{1,2}}, {{0,3}}, {{2,3}} }}, // 84
    {{ {{0,2}}, {{1,2}}, {{1,3}}, {{3,3}} }}, // 85
    {{ {{0,2}}, {{1,2}}, {{2,3}}, {{3,3}} }}, // 86
    {{ {{0,2}}, {{0,3}}, {{1,3}}, {{3,3}} }}, // 87
    {{ {{0,2}}, {{1,3}}, {{2,3}}, {{3,3}} }}, // 88
    {{ {{1,2}}, {{2,2}}, {{0,3}}, {{3,3}} }}, // 89
    {{ {{1,2}}, {{2,2}}, {{0,3}}, {{2,3}} }}, // 90
    {{ {{1,2}}, {{0,3}}, {{1,3}}, {{3,3}} }}, // 91
    {{ {{1,2}}, {{0,3}}, {{1,3}}, {{2,3}} }}, // 92
    {{ {{1,2}}, {{0,3}}, {{2,3}}, {{3,3}} }}  // 93
}};

constexpr std::array<std::array<unsigned char,8>,64> refinement_rules
{{
    {{  8, 99, 99, 99, 99, 99, 99, 99 }}, //  0
    {{ 21, 52, 99, 99, 99, 99, 99, 99 }}, //  1
    {{ 12, 37, 99, 99, 99, 99, 99, 99 }}, //  2
    {{ 24, 52, 62, 99, 99, 99, 99, 99 }}, //  3
    {{ 15, 27, 99, 99, 99, 99, 99, 99 }}, //  4
    {{ 21, 55, 71, 99, 99, 99, 99, 99 }}, //  5
    {{ 12, 40, 83, 99, 99, 99, 99, 99 }}, //  6
    {{ 24, 55, 65, 83, 99, 99, 99, 99 }}, //  7
    {{  9, 34, 99, 99, 99, 99, 99, 99 }}, //  8
    {{ 22, 52, 57, 99, 99, 99, 99, 99 }}, //  9
    {{ 13, 37, 43, 99, 99, 99, 99, 99 }}, // 10
    {{  0, 52, 62, 67, 99, 99, 99, 99 }}, // 11
    {{ 16, 28, 47, 89, 99, 99, 99, 99 }}, // 12
    {{ 22, 55, 57, 71, 99, 99, 99, 99 }}, // 13
    {{ 13, 40, 43, 83, 99, 99, 99, 99 }}, // 14
    {{  0, 55, 65, 67, 83, 99, 99, 99 }}, // 15
    {{ 10, 19, 99, 99, 99, 99, 99, 99 }}, // 16
    {{ 21, 53, 60, 99, 99, 99, 99, 99 }}, // 17
    {{ 14, 26, 38, 81, 99, 99, 99, 99 }}, // 18
    {{ 24, 53, 60, 62, 99, 99, 99, 99 }}, // 19
    {{ 17, 27, 30, 99, 99, 99, 99, 99 }}, // 20
    {{  1, 21, 71, 75, 99, 99, 99, 99 }}, // 21
    {{ 14, 26, 42, 83, 85, 99, 99, 99 }}, // 22
    {{  1, 24, 65, 75, 83, 99, 99, 99 }}, // 23
    {{  9, 35, 50, 99, 99, 99, 99, 99 }}, // 24
    {{ 22, 50, 53, 58, 99, 99, 99, 99 }}, // 25
    {{ 13, 38, 44, 81, 87, 99, 99, 99 }}, // 26
    {{  0,  3, 53, 63, 81, 87, 99, 99 }}, // 27
    {{ 16, 28, 48, 89, 91, 99, 99, 99 }}, // 28
    {{  1, 22, 72, 74, 89, 91, 99, 99 }}, // 29
    {{ 13, 42, 44, 83, 85, 87, 99, 99 }}, // 30
    {{  0,  1,  2,  3, 83, 85, 87, 99 }}, // 31
    {{ 11, 18, 99, 99, 99, 99, 99, 99 }}, // 32
    {{ 23, 25, 54, 56, 99, 99, 99, 99 }}, // 33
    {{ 12, 39, 46, 99, 99, 99, 99, 99 }}, // 34
    {{ 24, 54, 56, 64, 70, 99, 99, 99 }}, // 35
    {{ 15, 29, 32, 99, 99, 99, 99, 99 }}, // 36
    {{ 23, 25, 55, 73, 77, 99, 99, 99 }}, // 37
    {{  4, 12, 40, 86, 99, 99, 99, 99 }}, // 38
    {{  4, 24, 55, 65, 86, 99, 99, 99 }}, // 39
    {{  9, 36, 49, 99, 99, 99, 99, 99 }}, // 40
    {{ 22, 54, 56, 59, 79, 99, 99, 99 }}, // 41
    {{ 13, 39, 45, 49, 99, 99, 99, 99 }}, // 42
    {{  0, 54, 56, 64, 68, 79, 99, 99 }}, // 43
    {{ 16, 28, 47, 90, 93, 99, 99, 99 }}, // 44
    {{ 22, 55, 59, 73, 77, 79, 99, 99 }}, // 45
    {{  4, 13, 41, 47, 84, 93, 99, 99 }}, // 46
    {{  0,  4, 55, 66, 68, 77, 79, 99 }}, // 47
    {{ 10, 20, 33, 99, 99, 99, 99, 99 }}, // 48
    {{ 23, 25, 53, 61, 80, 99, 99, 99 }}, // 49
    {{ 14, 26, 38, 82, 88, 99, 99, 99 }}, // 50
    {{ 24, 53, 63, 69, 82, 88, 99, 99 }}, // 51
    {{ 17, 29, 31, 33, 99, 99, 99, 99 }}, // 52
    {{  1, 23, 25, 73, 76, 80, 99, 99 }}, // 53
    {{  4,  5, 14, 26, 42, 88, 99, 99 }}, // 54
    {{  1,  2,  4,  5, 24, 69, 88, 99 }}, // 55
    {{  7,  9, 35, 51, 99, 99, 99, 99 }}, // 56
    {{  7, 22, 51, 53, 58, 99, 99, 99 }}, // 57
    {{  6,  7, 13, 38, 44, 82, 99, 99 }}, // 58
    {{  0,  3,  6,  7, 53, 63, 82, 99 }}, // 59
    {{  7, 16, 28, 48, 90, 92, 99, 99 }}, // 60
    {{  1,  7, 22, 59, 73, 76, 78, 99 }}, // 61
    {{  4,  5,  6,  7, 13, 42, 44, 99 }}, // 62
    {{  0,  1,  2,  3,  4,  5,  6,  7 }}  // 63
}};

constexpr unsigned char refinement_children[64] = 
{
    0, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 4, 4, 4, 5, 2, 3, 4, 4, 3, 4, 5, 5, 3,
    4, 5, 6, 5, 6, 6, 7, 2, 4, 3, 5, 3, 5, 4, 5, 3, 5, 4, 6, 5, 6, 6, 7, 3, 5,
    5, 6, 4, 6, 6, 7, 4, 5, 6, 7, 6, 7, 7, 8
};

}

#endif

