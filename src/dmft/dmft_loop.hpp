#pragma once

#include "impurity_fine.hpp"
#include "impurity_raw.hpp"
#include "matsubara.hpp"
#include "parameters.hpp"
#include "random.hpp"

// Replaces legacy solver.cpp's bath_f2r/bath_r2f free functions and
// `struct DMFT_AF_solver`: the outer self-consistency loop coupling a
// RawImpuritySolver and a FineImpuritySolver on a two-sublattice
// antiferromagnetic square lattice (see lattice.hpp's dispersion()).
namespace dmft {

// Fits the raw solver's bath to the fine solver's current hybridization
// (Delta_u/Delta_d) at (fine.mu_fine, fine.hz_fine - lambda0*<Sz>), then
// diagonalizes it. Returns ln(Z). Equivalent to the original
// bath_f2r(solver_fine&, solver_raw&, double).
double bath_fine_to_raw(const FineImpuritySolver& fine, RawImpuritySolver& raw, Rng& rng,
                         double lambda0 = 0.0);

// Copies the raw solver's converged state into the fine solver (calc_g so
// its derivatives are refreshed) and immediately expands g0/gz/gu/gd/...
// around it via set_g(). Equivalent to the original
// bath_r2f(solver_raw&, solver_fine&). (Unused by AntiferromagneticDMFTLoop
// below -- same as in the original -- ported for interface parity.)
void bath_raw_to_fine(const RawImpuritySolver& raw, FineImpuritySolver& fine, Rng& rng);

class AntiferromagneticDMFTLoop {
public:
    AntiferromagneticDMFTLoop(const Parameters& params, const MatsubaraGrid& grid,
                               double lambda0 = 0.0);

    // Seeds Sf.Delta_u with a single-pole guess (see impurity_fine.cpp's
    // note on the original's Delta_d never being touched here -- preserved
    // as in the original, not fixed), sets Sf.mu_fine/hz_fine, and fits Sr
    // to it. Equivalent to the original init(double hz=0).
    void init(double hz, Rng& rng);

    // One step of the lattice self-consistency condition (a two-sublattice
    // Dyson equation, folding the reduced Brillouin zone via kx+Lx/2,
    // ky+Ly/2), updating Sf.Delta_u/Delta_d in place by `factor` of the
    // mismatch. Equivalent to the original iter_Delta(double).
    void iter_Delta(double factor);

    // Iterates bath_fine_to_raw() + Sf.init(Sr, calc_g=false) +
    // iter_Delta() up to 1000 times, stopping once Sf.gu[0] stops moving
    // by more than 1e-4. Returns the stopping iteration, or -1 if it
    // never converged. Equivalent to the original raw_loops(double).
    int raw_loops(double factor, Rng& rng);

    // Iterates Sf.set_g() + iter_Delta() + an hz_fine Lagrange-multiplier
    // update (for the lambda0 staggered-field constraint) up to 1000
    // times, stopping once Sf.gu[0] stops moving by more than 1e-6.
    // Returns the stopping iteration, or -1 if it never converged.
    // Equivalent to the original fine_loops(double).
    int fine_loops(double factor, Rng& rng);

    RawImpuritySolver Sr;
    FineImpuritySolver Sf;
    bool af_flag = true;
    double lambda0() const { return lambda0_; }

private:
    const Parameters* params_;
    const MatsubaraGrid* grid_;
    double lambda0_;
};

}  // namespace dmft
