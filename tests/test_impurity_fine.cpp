// Cross-checks src/dmft/impurity_fine.* (FineImpuritySolver) against the
// original solver_fine (via tests/impurity_fine_bridge.cpp). Builds a raw
// solver's converged state once (already validated bit-exact against
// legacy in test_impurity_raw.cpp), then compares FineImpuritySolver's
// deterministic set_g_raw() and its RNG-driven set_g() -- called several
// times in a row, exercising expand()'s lazy-once basis computation and
// its reuse on subsequent calls -- against the legacy solver_fine.
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

#include "../src/dmft/impurity_fine.hpp"
#include "../src/dmft/impurity_raw.hpp"
#include "impurity_fine_bridge.hpp"

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

void compare_snapshot(const char* prefix, const dmft::FineImpuritySolver& fine,
                       const legacy_bridge::FineSnapshot& legacy, double tol) {
    char label[160];
    struct Field {
        const char* name;
        const std::vector<Complex>* new_field;
        const std::vector<Complex>* legacy_field;
    };
    const Field fields[] = {
        {"g0", &fine.g0, &legacy.g0}, {"gz", &fine.gz, &legacy.gz},
        {"gu", &fine.gu, &legacy.gu}, {"gd", &fine.gd, &legacy.gd},
        {"gx", &fine.gx, &legacy.gx}, {"gy", &fine.gy, &legacy.gy},
        {"gax", &fine.gax, &legacy.gax}, {"gay", &fine.gay, &legacy.gay},
    };
    for (const auto& f : fields) {
        std::snprintf(label, sizeof(label), "%s: %s matches legacy solver_fine", prefix, f.name);
        expect(max_abs_diff(*f.new_field, *f.legacy_field) < tol, label);
    }
}

}  // namespace

int main() {
    dmft::Parameters params;  // defaults match legacy/parameters.cpp exactly
    dmft::MatsubaraGrid grid(params.Nw, params.beta);

    // Same well-conditioned 3-pole hybridization as test_impurity_raw.cpp.
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
    const double mu_loc = 0.1, h_loc = 0.0;

    const unsigned int seed_raw = 20260903;
    legacy_bridge::reset_solvers();
    legacy_bridge::reseed(seed_raw);
    legacy_bridge::raw_then_fine_init(delta_up, delta_down, mu_loc, h_loc, /*calc_g=*/true);

    dmft::RawImpuritySolver raw(params, grid);
    dmft::Rng rng_raw(seed_raw);
    raw.init(delta_up, delta_down, mu_loc, h_loc, rng_raw);
    dmft::FineImpuritySolver fine(params, grid);
    fine.init(raw, /*calc_g=*/true);

    const double tol = 1e-6;

    {
        const legacy_bridge::FineSnapshot legacy = legacy_bridge::fine_set_g_raw();
        fine.set_g_raw();
        compare_snapshot("set_g_raw", fine, legacy, tol);
    }

    // Small perturbations away from the raw-solver state, well inside the
    // linear-response regime get_derivatives() was built for.
    std::vector<Complex> delta_u(params.Nw), delta_d(params.Nw), delta_x(params.Nw),
        delta_y(params.Nw), delta_ax(params.Nw), delta_ay(params.Nw);
    for (int w = 0; w < params.Nw; ++w) {
        delta_u[w] = 0.01 / (grid.iw(w) - 0.15);
        delta_d[w] = 0.008 / (grid.iw(w) + 0.1);
        delta_x[w] = 0.002 / (grid.iw(w) - 0.05);
        delta_y[w] = 0.0015 / (grid.iw(w) + 0.05);
        delta_ax[w] = 0.001 / (grid.iw(w) - 0.02);
        delta_ay[w] = 0.0012 / (grid.iw(w) + 0.02);
    }
    const double mu_fine = mu_loc + 0.005, hz_fine = 0.003, hx_fine = 0.002, hax_fine = 0.001,
                 hay_fine = 0.0015;

    legacy_bridge::fine_set_deltas(delta_u, delta_d, delta_x, delta_y, delta_ax, delta_ay, mu_fine,
                                    hz_fine, hx_fine, hax_fine, hay_fine);
    fine.Delta_u = delta_u;
    fine.Delta_d = delta_d;
    fine.Delta_x = delta_x;
    fine.Delta_y = delta_y;
    fine.Delta_ax = delta_ax;
    fine.Delta_ay = delta_ay;
    fine.mu_fine = mu_fine;
    fine.hz_fine = hz_fine;
    fine.hx_fine = hx_fine;
    fine.hax_fine = hax_fine;
    fine.hay_fine = hay_fine;

    const unsigned int seed_fine = 777;
    legacy_bridge::reseed(seed_fine);
    dmft::Rng rng_fine(seed_fine);

    const int n_calls = 3;  // first call triggers expand()'s lazy basis setup
    for (int call = 0; call < n_calls; ++call) {
        const legacy_bridge::FineSnapshot legacy = legacy_bridge::fine_set_g();
        fine.set_g(rng_fine);
        char prefix[32];
        std::snprintf(prefix, sizeof(prefix), "set_g call %d", call);
        compare_snapshot(prefix, fine, legacy, tol);
    }

    std::printf("\n%s\n", g_failures == 0 ? "ALL TESTS PASSED" : "TESTS FAILED");
    return g_failures == 0 ? 0 : 1;
}
