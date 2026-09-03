#pragma once

// Plain, macro-free bridge into the legacy solver_raw + solver_fine pair,
// for cross-checking src/dmft/impurity_fine.* against the original. See
// tests/legacy_bridge.hpp for why this split exists.

#include <complex>
#include <vector>

namespace legacy_bridge {

using Complex = std::complex<double>;

// Reseeds the legacy global RNG (INT_RANDOM). Call once before a sequence
// of calls meant to be compared against a freshly Rng-seeded sequence on
// the new-code side.
void reseed(unsigned int seed);

// (Re)constructs the one legacy solver_raw + solver_fine instance pair
// that the functions below operate on.
void reset_solvers();

// Sr.init(delta_up, delta_down, mu_loc, h_loc), then Sf.init(Sr, calc_g).
void raw_then_fine_init(const std::vector<Complex>& delta_up, const std::vector<Complex>& delta_down,
                         double mu_loc, double h_loc, bool calc_g);

struct FineSnapshot {
    std::vector<Complex> g0, gz, gu, gd, gx, gy, gax, gay;
};

// Sf.set_g_raw(), then a snapshot of Sf.g0/gz/gu/gd/gx/gy/gax/gay.
FineSnapshot fine_set_g_raw();

// Sets Sf.Delta_u/Delta_d/Delta_x/Delta_y/Delta_ax/Delta_ay and
// Sf.mu_fine/hz_fine/hx_fine/hax_fine/hay_fine (Sf.hy_fine is left
// untouched: the original's set_g() never reads it -- gy is computed with
// hx_fine, same as gx; see impurity_fine.cpp's set_g() for where this is
// preserved rather than "fixed").
void fine_set_deltas(const std::vector<Complex>& delta_u, const std::vector<Complex>& delta_d,
                      const std::vector<Complex>& delta_x, const std::vector<Complex>& delta_y,
                      const std::vector<Complex>& delta_ax, const std::vector<Complex>& delta_ay,
                      double mu_fine, double hz_fine, double hx_fine, double hax_fine,
                      double hay_fine);

// Sf.set_g(), then a snapshot of Sf.g0/gz/gu/gd/gx/gy/gax/gay.
FineSnapshot fine_set_g();

}  // namespace legacy_bridge
