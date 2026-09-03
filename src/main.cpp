#include <iostream>

#include "dmft/dmft_loop.hpp"

// Replaces legacy/main.cpp: runs the same two-sublattice antiferromagnetic
// DMFT calculation (DMFT_AF_solver(0.); solver.init(.003);
// solver.raw_loops(1); solver.fine_loops(1);) against the rewritten
// dmft:: classes, and prints the same diagnostics.
int main() {
    dmft::Parameters params;
    dmft::MatsubaraGrid grid(params.Nw, params.beta);
    dmft::AntiferromagneticDMFTLoop loop(params, grid, /*lambda0=*/0.0);
    dmft::Rng rng;  // time-seeded, matching the original's default generator

    loop.init(0.003, rng);
    const double field_u = loop.Sf.hz_fine - loop.Sf.magnetization() * loop.lambda0();
    std::cout << loop.raw_loops(1.0, rng) << "  "
              << (1.0 / params.beta) * loop.Sf.magnetization() / field_u << "  s=" << loop.Sf.magnetization()
              << "  at h=" << field_u << "\n";

    std::cout << loop.fine_loops(1.0, rng) << "  "
              << (1.0 / params.beta) * loop.Sf.magnetization() /
                     (loop.Sf.hz_fine - loop.Sf.magnetization() * loop.lambda0())
              << "\n";

    return 0;
}
