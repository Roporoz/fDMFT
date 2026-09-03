// Cross-checks src/dmft/impurity_raw.* (RawImpuritySolver) against the
// original solver_raw (via tests/impurity_raw_bridge.cpp) across several
// successive init() calls on one instance of each -- exercising the same
// RNG-driven bath search, and the same up-then-down adaptive state
// (trial budget / restart jump size) coupling between calls, that the real
// DMFT driver relies on.
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

#include "../src/dmft/impurity_raw.hpp"
#include "impurity_raw_bridge.hpp"

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

double max_abs_diff(const std::vector<double>& a, const std::vector<double>& b) {
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
    return m;
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

    // A 3-pole hybridization (matching Nbath_raw=3) rather than a single
    // pole: fitting one pole with 3 free bath levels is a degenerate,
    // ill-conditioned problem (many equally-good, some unphysical, local
    // minima), which makes the search burn through its whole trial budget
    // on essentially every call -- true of the original algorithm just as
    // much as this rewrite, but a poor choice for a fast, well-posed test.
    std::vector<Complex> delta_up(params.Nw), delta_down(params.Nw);
    const double eu_poles[3] = {-0.5, 0.1, 0.6};
    const double t2u_poles[3] = {0.05, 0.08, 0.03};
    const double ed_poles[3] = {-0.4, 0.05, 0.55};
    const double t2d_poles[3] = {0.04, 0.07, 0.035};
    for (int w = 0; w < params.Nw; ++w) {
        delta_up[w] = 0.0;
        delta_down[w] = 0.0;
        for (int l = 0; l < 3; ++l) {
            delta_up[w] += t2u_poles[l] / (grid.iw(w) - eu_poles[l]);
            delta_down[w] += t2d_poles[l] / (grid.iw(w) - ed_poles[l]);
        }
    }

    const unsigned int seed = 20260903;
    legacy_bridge::reset_raw_solver(seed);
    dmft::RawImpuritySolver new_solver(params, grid);
    dmft::Rng rng(seed);

    const int n_iterations = 5;
    const double tol = 1e-6;
    for (int iter = 0; iter < n_iterations; ++iter) {
        const double mu_loc = 0.1;
        const double h_loc = 0.01 * iter;

        const legacy_bridge::RawSolverSnapshot legacy =
            legacy_bridge::raw_solver_init(delta_up, delta_down, mu_loc, h_loc);
        const double ln_z_new = new_solver.init(delta_up, delta_down, mu_loc, h_loc, rng);

        char label[128];
        std::snprintf(label, sizeof(label), "iter %d: ln(Z) matches legacy solver_raw::init()", iter);
        expect(std::abs(ln_z_new - legacy.ln_z) < tol, label);

        std::snprintf(label, sizeof(label), "iter %d: tolerance (bath residual) matches", iter);
        expect(std::abs(new_solver.tolerance - legacy.tolerance) < tol, label);

        std::snprintf(label, sizeof(label), "iter %d: t2u/eu (up bath) match", iter);
        expect(max_abs_diff(new_solver.t2u, legacy.t2u) < tol &&
                   max_abs_diff(new_solver.eu, legacy.eu) < tol,
               label);

        std::snprintf(label, sizeof(label), "iter %d: t2d/ed (down bath) match", iter);
        expect(max_abs_diff(new_solver.t2d, legacy.t2d) < tol &&
                   max_abs_diff(new_solver.ed, legacy.ed) < tol,
               label);

        std::snprintf(label, sizeof(label), "iter %d: gu/gd (impurity Green's function) match", iter);
        expect(max_abs_diff(new_solver.gu, legacy.gu) < tol &&
                   max_abs_diff(new_solver.gd, legacy.gd) < tol,
               label);

        std::snprintf(label, sizeof(label), "iter %d: sigmau/sigmad (self-energy) match", iter);
        expect(max_abs_diff(new_solver.sigmau, legacy.sigmau) < tol &&
                   max_abs_diff(new_solver.sigmad, legacy.sigmad) < tol,
               label);
    }

    std::printf("\n%s\n", g_failures == 0 ? "ALL TESTS PASSED" : "TESTS FAILED");
    return g_failures == 0 ? 0 : 1;
}
