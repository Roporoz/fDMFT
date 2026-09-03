#pragma once

// Plain, macro-free bridge into the legacy solver_raw class, for
// cross-checking src/dmft/impurity_raw.* against it. See
// tests/legacy_bridge.hpp for why this split exists (this header is
// included from test code built under C++17; the .cpp implementing it is
// the only place that pulls in legacy/headers.new.h + parameters.cpp +
// solver.cpp, under -std=gnu++14, exactly as legacy/main.cpp does).

#include <complex>
#include <vector>

namespace legacy_bridge {

using Complex = std::complex<double>;

struct RawSolverSnapshot {
    double ln_z = 0.0;
    double tolerance = 0.0;
    std::vector<double> t2u, t2d, eu, ed;
    std::vector<Complex> gu, gd;
    std::vector<Complex> sigmau, sigmad;
};

// (Re)constructs the one legacy solver_raw instance that raw_solver_init()
// calls into, and reseeds INT_RANDOM=seed. Call once before a sequence of
// raw_solver_init() calls meant to be compared against a freshly
// constructed dmft::RawImpuritySolver seeded the same way.
void reset_raw_solver(unsigned int seed);

// Calls solver_raw::init(delta_up, delta_down, mu_loc, h_loc) on the
// instance from the last reset_raw_solver(), returning a snapshot of its
// public state afterwards. delta_up/delta_down must have Nw entries
// (legacy/parameters.cpp's Nw=200).
RawSolverSnapshot raw_solver_init(const std::vector<Complex>& delta_up,
                                   const std::vector<Complex>& delta_down, double mu_loc,
                                   double h_loc);

}  // namespace legacy_bridge
