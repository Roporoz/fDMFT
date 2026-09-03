#pragma once

#include <cmath>

#include "core.hpp"
#include "parameters.hpp"

// Replaces legacy solver.cpp's free function eps(int, int): the tight-
// binding dispersion of a square lattice with nearest-neighbor hopping t
// and two diagonal next-nearest-neighbor hoppings t_1 (along kx+ky) and
// t_2 (along kx-ky) -- e.g. for a bilayer or next-nearest-neighbor model.
namespace dmft {

inline double dispersion(const Parameters& params, int kx, int ky) {
    const double qx = kx * 2.0 * pi / params.Lx;
    const double qy = ky * 2.0 * pi / params.Ly;
    return -2.0 * params.t * (std::cos(qx) + std::cos(qy)) - 2.0 * params.t_1 * std::cos(qx + qy) -
           2.0 * params.t_2 * std::cos(qx - qy);
}

}  // namespace dmft
