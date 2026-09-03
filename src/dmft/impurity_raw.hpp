#pragma once

#include <vector>

#include "core.hpp"
#include "linalg.hpp"
#include "matsubara.hpp"
#include "parameters.hpp"
#include "random.hpp"

// Replaces legacy/solver.cpp's `class solver_raw`: an exact-diagonalization
// solver for a two-orbital (spin up/down) single-impurity Anderson model
// with a finite discrete bath, fit by chi^2 minimization against a given
// hybridization function.
//
// Deliberately NOT ported here (out of scope for this pass): set_g2() and
// the g2tau_*/gamma_* two-particle vertex machinery. main.cpp's actual
// DMFT driver never calls them, so there is no exercised behavior to
// validate the port against yet; they're tracked as follow-up work.
namespace dmft {

class RawImpuritySolver {
public:
    RawImpuritySolver(const Parameters& params, const MatsubaraGrid& grid);

    // Fits both spin baths to (delta_up, delta_down) -- see fit_bath() --
    // then diagonalizes the resulting impurity model at (mu_loc, h_loc) --
    // see diagonalize(). Returns ln(Z) (matching the original's return
    // value); the achieved bath-fit residual is left in `tolerance`.
    // Equivalent to the original solver_raw::init().
    double init(const std::vector<Complex>& delta_up, const std::vector<Complex>& delta_down,
                double mu_loc, double h_loc, Rng& rng);

    // Random-search + Newton-refinement fit of both spin baths (t2u/eu and
    // t2d/ed) to (delta_up, delta_down), by chi^2 minimization. Returns the
    // summed residual (what `tolerance` holds after init()).
    double fit_bath(const std::vector<Complex>& delta_up, const std::vector<Complex>& delta_down,
                     Rng& rng);

    // Builds the impurity Hamiltonian from the *current* t2u/t2d/eu/ed at
    // chemical potential mu_loc and field h_loc, exactly diagonalizes it,
    // and computes gu/gd/sigmau/sigmad from the resulting eigenbasis.
    // Returns ln(Z). Equivalent to the tail of the original init(), after
    // its two adjust_bath_levels() calls.
    double diagonalize(double mu_loc, double h_loc);

    Complex Gu(int w) const { return w >= 0 ? gu[w] : std::conj(gu[-w - 1]); }
    Complex Gd(int w) const { return w >= 0 ? gd[w] : std::conj(gd[-w - 1]); }

    // Equivalent to the original s(): 2*<Sz>.
    double magnetization() const;
    // Equivalent to the original n(): <n_up> + <n_down>.
    double occupation() const;

    Complex Gudh(int w) const { return w >= 0 ? gudh[w] : std::conj(gudh[-w - 1]); }
    Complex Gddh(int w) const { return w >= 0 ? gddh[w] : std::conj(gddh[-w - 1]); }
    Complex Gxdh(int w) const { return w >= 0 ? gxdh[w] : std::conj(gxdh[-w - 1]); }
    // Equivalent to the original dszdh()/dsxdh(). NOTE: in the original,
    // gudh/gddh/gxdh are allocated but never written to by anything --
    // solver_raw::set_g2() (not ported; see above) is the only method that
    // could plausibly fill them, and even it doesn't. These two accessors
    // and the three above are ported for interface parity but currently
    // always read back zero.
    double dm_dh_z() const;
    double dm_dh_x() const;

    // Reconstructs the hybridization function from the fitted bath
    // (delta_u_out/delta_d_out), alongside a copy of gu/gd. Equivalent to
    // the original store(). (Unused by the current DMFT driver, ported for
    // interface parity.)
    void store(std::vector<Complex>& delta_u_out, std::vector<Complex>& delta_d_out,
               std::vector<Complex>& g_u_out, std::vector<Complex>& g_d_out) const;

    double tolerance = 0.0;
    std::vector<double> t2u, t2d, eu, ed;
    double mu_loc = 0.0, h_loc = 0.0;
    std::vector<Complex> gu, gd;
    std::vector<Complex> sigmau, sigmad;
    std::vector<Complex> gudh, gddh, gxdh;

private:
    static int occ(int state, int site) { return (state >> site) & 1; }

    // Given fixed bath levels e_bath, solves in closed form for the
    // couplings t2_bath minimizing the chi^2 residual between delta_w and
    // sum_i t2_bath[i]/(iw - e_bath[i]), and returns that residual. Pure
    // (same e_bath and delta_w always give the same answer) and RNG-free.
    // Equivalent to the original adjust_bath_couplings().
    double adjust_bath_couplings(const std::vector<Complex>& delta_w, std::vector<double>& t2_bath,
                                  const std::vector<double>& e_bath) const;

    // Random-restart, Newton-refined search over bath levels e_bath for
    // the fit minimizing adjust_bath_couplings()'s residual. Adapts its
    // own trial budget and restart jump size across successive calls on
    // this solver (trial_budget_/jump_scale_ below), the same way the
    // original's function-local statics did -- except scoped to this
    // instance rather than shared, unscoped, across every solver_raw in
    // the process (harmless in practice, since the actual driver only
    // ever constructs one). Equivalent to adjust_bath_levels().
    double adjust_bath_levels(const std::vector<Complex>& delta_w, std::vector<double>& t2_bath,
                               std::vector<double>& e_bath, Rng& rng);

    const Parameters* params_;
    const MatsubaraGrid* grid_;
    int n_fock_;

    Matrix<double> H_, Psi_, cu_psi_, cd_psi_, cu_, cd_;
    std::vector<double> E_;
    std::vector<int> nud_;
    long double Z_ = 0.0;

    int trial_budget_ = 10000;
    double jump_scale_ = 1.0;

    // Reused scratch space for adjust_bath_couplings()/adjust_bath_levels():
    // the former is called millions of times over a full bath search (up to
    // trial_budget_ trials x up to 50 Newton steps x ~2*Nbath_raw+3 calls
    // each), so allocating a fresh Nbath_raw x Nbath_raw matrix per call --
    // as opposed to the original's single reused `static` buffer -- is a
    // real, measured slowdown, not just style. Mutable because
    // adjust_bath_couplings() is logically const (same inputs always give
    // the same answer); this is scratch space, not observable state.
    mutable Matrix<double> bath_fit_matrix_;
    Matrix<double> newton_hessian_;
    std::vector<double> newton_gradient_;
};

}  // namespace dmft
