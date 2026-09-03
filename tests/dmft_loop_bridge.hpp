#pragma once

// Plain, macro-free bridge into the legacy DMFT_AF_solver, for
// cross-checking src/dmft/dmft_loop.* against the original. See
// tests/legacy_bridge.hpp for why this split exists.

#include <complex>
#include <vector>

namespace legacy_bridge {

using Complex = std::complex<double>;

struct DmftLoopSnapshot {
    std::vector<Complex> gu, gd, delta_u, delta_d;
    double hz_fine = 0.0;
    int raw_loops_result = 0;
    int fine_loops_result = 0;
};

// Reseeds INT_RANDOM=seed, then constructs DMFT_AF_solver(lambda0) and
// runs exactly what legacy/main.cpp runs: .init(hz), .raw_loops(factor),
// .fine_loops(factor). Returns a snapshot of the final state.
DmftLoopSnapshot run(double lambda0, double hz, double factor, unsigned int seed);

}  // namespace legacy_bridge
