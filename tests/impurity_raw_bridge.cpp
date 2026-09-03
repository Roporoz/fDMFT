#include "impurity_raw_bridge.hpp"

#include "../legacy/headers.new.h"
#include "../legacy/parameters.cpp"
#include "../legacy/solver.cpp"

namespace legacy_bridge {

namespace {
solver_raw* g_solver = nullptr;
}

void reset_raw_solver(unsigned int seed) {
    delete g_solver;
    g_solver = new solver_raw();
    INT_RANDOM = static_cast<int>(seed);
}

RawSolverSnapshot raw_solver_init(const std::vector<Complex>& delta_up,
                                   const std::vector<Complex>& delta_down, double mu_loc,
                                   double h_loc) {
    complex* du = new complex[Nw];
    complex* dd = new complex[Nw];
    for (int w = 0; w < Nw; ++w) {
        du[w] = delta_up[w];
        dd[w] = delta_down[w];
    }

    RawSolverSnapshot snap;
    snap.ln_z = g_solver->init(du, dd, mu_loc, h_loc);
    snap.tolerance = g_solver->tolerance;
    snap.t2u.assign(g_solver->t2u, g_solver->t2u + Nbath_raw);
    snap.t2d.assign(g_solver->t2d, g_solver->t2d + Nbath_raw);
    snap.eu.assign(g_solver->eu, g_solver->eu + Nbath_raw);
    snap.ed.assign(g_solver->ed, g_solver->ed + Nbath_raw);
    snap.gu.assign(g_solver->gu, g_solver->gu + Nw);
    snap.gd.assign(g_solver->gd, g_solver->gd + Nw);
    snap.sigmau.assign(g_solver->sigmau, g_solver->sigmau + Nw);
    snap.sigmad.assign(g_solver->sigmad, g_solver->sigmad + Nw);

    delete[] du;
    delete[] dd;
    return snap;
}

}  // namespace legacy_bridge
