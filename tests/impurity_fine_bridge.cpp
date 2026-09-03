#include "impurity_fine_bridge.hpp"

#include "../legacy/headers.new.h"
#include "../legacy/parameters.cpp"
#include "../legacy/solver.cpp"

namespace legacy_bridge {

namespace {
solver_raw* g_raw = nullptr;
solver_fine* g_fine = nullptr;

FineSnapshot snapshot() {
    FineSnapshot snap;
    snap.g0.assign(g_fine->g0, g_fine->g0 + Nw);
    snap.gz.assign(g_fine->gz, g_fine->gz + Nw);
    snap.gu.assign(g_fine->gu, g_fine->gu + Nw);
    snap.gd.assign(g_fine->gd, g_fine->gd + Nw);
    snap.gx.assign(g_fine->gx, g_fine->gx + Nw);
    snap.gy.assign(g_fine->gy, g_fine->gy + Nw);
    snap.gax.assign(g_fine->gax, g_fine->gax + Nw);
    snap.gay.assign(g_fine->gay, g_fine->gay + Nw);
    return snap;
}
}  // namespace

void reseed(unsigned int seed) { INT_RANDOM = static_cast<int>(seed); }

void reset_solvers() {
    delete g_raw;
    delete g_fine;
    g_raw = new solver_raw();
    g_fine = new solver_fine();
}

void raw_then_fine_init(const std::vector<Complex>& delta_up, const std::vector<Complex>& delta_down,
                         double mu_loc, double h_loc, bool calc_g) {
    complex* du = new complex[Nw];
    complex* dd = new complex[Nw];
    for (int w = 0; w < Nw; ++w) {
        du[w] = delta_up[w];
        dd[w] = delta_down[w];
    }
    g_raw->init(du, dd, mu_loc, h_loc);
    delete[] du;
    delete[] dd;

    g_fine->init(*g_raw, calc_g);
}

FineSnapshot fine_set_g_raw() {
    g_fine->set_g_raw();
    return snapshot();
}

void fine_set_deltas(const std::vector<Complex>& delta_u, const std::vector<Complex>& delta_d,
                      const std::vector<Complex>& delta_x, const std::vector<Complex>& delta_y,
                      const std::vector<Complex>& delta_ax, const std::vector<Complex>& delta_ay,
                      double mu_fine, double hz_fine, double hx_fine, double hax_fine,
                      double hay_fine) {
    for (int w = 0; w < Nw; ++w) {
        g_fine->Delta_u[w] = delta_u[w];
        g_fine->Delta_d[w] = delta_d[w];
        g_fine->Delta_x[w] = delta_x[w];
        g_fine->Delta_y[w] = delta_y[w];
        g_fine->Delta_ax[w] = delta_ax[w];
        g_fine->Delta_ay[w] = delta_ay[w];
    }
    g_fine->mu_fine = mu_fine;
    g_fine->hz_fine = hz_fine;
    g_fine->hx_fine = hx_fine;
    g_fine->hax_fine = hax_fine;
    g_fine->hay_fine = hay_fine;
}

FineSnapshot fine_set_g() {
    g_fine->set_g();
    return snapshot();
}

}  // namespace legacy_bridge
