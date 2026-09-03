// Cross-checks src/dmft/dmft_loop.* (AntiferromagneticDMFTLoop) against
// the original DMFT_AF_solver end to end -- init(), raw_loops(), and
// fine_loops(), exactly as legacy/main.cpp runs them -- rather than
// against any already-validated piece in isolation. This is the top-level
// integration test for the whole rewrite: RawImpuritySolver and
// FineImpuritySolver were each already validated bit-exact on their own
// (see test_impurity_raw.cpp/test_impurity_fine.cpp); what's new and
// exercised only here is the coupling logic itself -- iter_Delta()'s
// lattice self-consistency step and bath_fine_to_raw()'s argument wiring.
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

#include "../src/dmft/dmft_loop.hpp"
#include "dmft_loop_bridge.hpp"

namespace {

using dmft::Complex;

int g_failures = 0;

void expect(bool condition, const char* what) {
    if (!condition) {
        std::printf("FAIL: %s\n", what);
        ++g_failures;
    } else {
        std::printf("ok:   %s\n", what);
    }
}

double max_abs_diff(const std::vector<Complex>& a, const std::vector<Complex>& b) {
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
    return m;
}

}  // namespace

int main() {
    dmft::Parameters params;  // defaults match legacy/parameters.cpp exactly
    dmft::MatsubaraGrid grid(params.Nw, params.beta);

    // Exactly legacy/main.cpp's call: DMFT_AF_solver solver(0.);
    // solver.init(.003); solver.raw_loops(1); solver.fine_loops(1);
    const double lambda0 = 0.0, hz = 0.003, factor = 1.0;
    const unsigned int seed = 20260903;

    const legacy_bridge::DmftLoopSnapshot legacy = legacy_bridge::run(lambda0, hz, factor, seed);

    dmft::AntiferromagneticDMFTLoop loop(params, grid, lambda0);
    dmft::Rng rng(seed);
    loop.init(hz, rng);
    const int raw_result = loop.raw_loops(factor, rng);
    const int fine_result = loop.fine_loops(factor, rng);

    const double tol = 1e-6;

    std::printf("raw_loops: new=%d legacy=%d; fine_loops: new=%d legacy=%d\n", raw_result,
                legacy.raw_loops_result, fine_result, legacy.fine_loops_result);

    expect(raw_result == legacy.raw_loops_result, "raw_loops() converges at the same iteration");
    expect(fine_result == legacy.fine_loops_result, "fine_loops() converges at the same iteration");
    expect(std::abs(loop.Sf.hz_fine - legacy.hz_fine) < tol, "final hz_fine matches");
    expect(max_abs_diff(loop.Sf.gu, legacy.gu) < tol, "final Sf.gu matches legacy DMFT_AF_solver");
    expect(max_abs_diff(loop.Sf.gd, legacy.gd) < tol, "final Sf.gd matches legacy DMFT_AF_solver");
    expect(max_abs_diff(loop.Sf.Delta_u, legacy.delta_u) < tol,
           "final Sf.Delta_u matches legacy DMFT_AF_solver");
    expect(max_abs_diff(loop.Sf.Delta_d, legacy.delta_d) < tol,
           "final Sf.Delta_d matches legacy DMFT_AF_solver");

    std::printf("\n%s\n", g_failures == 0 ? "ALL TESTS PASSED" : "TESTS FAILED");
    return g_failures == 0 ? 0 : 1;
}
