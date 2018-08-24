#ifndef BEM_LATTICE_H
#define BEM_LATTICE_H

#include "types.h"
#include "bary.h"
#include "geometry/point.h"

#include <array>

namespace bem
{

namespace lattice
{

template <uint order>
constexpr uint size();

template <uint order>
using positions = std::array<bary,size<order>()>;

template <uint order>
using index = std::array<uint,size<order>()>;

template <uint order>
using nodes = std::array<geometry::point,size<order>()>;

template <uint order>
constexpr positions<order> lattice();

}

}

#include "lattice.tpp"
#endif

