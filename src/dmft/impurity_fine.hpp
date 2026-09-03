#pragma once

#include <vector>

#include "core.hpp"
#include "impurity_raw.hpp"
#include "linalg.hpp"
#include "matsubara.hpp"
#include "parameters.hpp"
#include "random.hpp"

// Replaces legacy/solver.cpp's `class solver_fine`: a perturbative
// (linear-response) correction to RawImpuritySolver's Green's functions,
// built from finite-difference derivatives of the same kind of exact
// diagonalization around the raw solver's converged state. Where
// RawImpuritySolver re-diagonalizes exactly on every DMFT iteration,
// FineImpuritySolver diagonalizes once (in init()) to get those
// derivatives, then reuses them cheaply across many set_g() calls with a
// changing hybridization/field (see DMFT_AF_solver::fine_loops in the
// original, ported next).
namespace dmft {

class FineImpuritySolver {
public:
    FineImpuritySolver(const Parameters& params, const MatsubaraGrid& grid);

    // Copies the raw solver's current bath (t2u/t2d/eu/ed, mu_loc, h_loc)
    // and one-particle Green's functions into g0/gz/gu/gd (gx/gy/gax/gay
    // reset to zero). If calc_g, also (re)computes the field/hybridization
    // derivatives that set_g() uses -- an expensive step (dozens of exact
    // diagonalizations; see get_derivatives()) meant to be done once per
    // raw-solver state, not on every set_g() call.
    // Equivalent to the original init(solver_raw&, bool).
    void init(const RawImpuritySolver& raw, bool calc_g = true);

    double magnetization() const;  // was sz()
    double occupation() const;     // was n()

    // Sets g0/gz/gu/gd to the raw solver's Green's functions directly (and
    // gx/gy/gax/gay to zero), ignoring Delta_u/Delta_d/... and mu_fine/
    // hz_fine/... entirely. Equivalent to the original set_g_raw().
    void set_g_raw();

    // Sets g0/gz/gu/gd/gx/gy/gax/gay via the linear-response expansion
    // around the raw solver's state (from the last init()), using the
    // *current* Delta_u/Delta_d/Delta_x/Delta_y/Delta_ax/Delta_ay and
    // mu_fine/hz_fine/hx_fine/hax_fine/hay_fine. Equivalent to the
    // original set_g(). Draws random restart vectors from `rng`, but only
    // on this solver's very first call ever (see expand()) -- pass the
    // same Rng across a sequence of set_g() calls the way the original
    // read from one shared global generator.
    void set_g(Rng& rng);

    std::vector<Complex> g0, gz, gu, gd, gx, gy, gax, gay;
    std::vector<Complex> Delta_u, Delta_d, Delta_x, Delta_y, Delta_ax, Delta_ay;
    double mu_fine = 0.0, hz_fine = 0.0, hx_fine = 0.0, hy_fine = 0.0, hax_fine = 0.0,
           hay_fine = 0.0;

private:
    static int occ(int state, int site) { return (state >> site) & 1; }

    // Diagonalizes the same-shaped impurity Hamiltonian as
    // RawImpuritySolver::diagonalize() (no quantum-number blocking here,
    // since the anomal=true perturbation below doesn't conserve total
    // particle number), with a field perturbation dh added either to
    // (n_up - n_down) [anomal=true] or (n_up + n_down) [anomal=false].
    // Fills gu_/gd_/gx_/ga_ from the resulting eigenbasis.
    // Equivalent to the original get_g_static(double, bool).
    void get_g_static(double dh, bool anomal);

    // As above, but on a Fock space with one extra bath level added at
    // energy e_probe, coupled by t_ to spin up (flag=0), spin down
    // (flag=1), or both with an anomalous sign (flag=2).
    // Equivalent to the original get_g(double, double, int).
    void get_g(double t_, double e_probe, int flag);

    // Finite-difference derivatives of gu_/gd_/gx_/ga_ at the raw state:
    // w.r.t. mu and hz (via get_g_static at +-dh/2), and w.r.t. the
    // hybridization function at each of the Nbath_fine probe energies e_
    // (via get_g at two step sizes, Richardson-extrapolated). Fills
    // dg0_mu_/dg0_hz_/.../dg0_D0_/.../dga_Da_.
    // Equivalent to the original get_derivatives().
    void get_derivatives();

    double scalar(const std::vector<Complex>& r1, const std::vector<Complex>& r2) const;

    // Computes, once per solver (lazily -- see basis_ready_ -- since it's
    // an expensive Eigen_leading call that doesn't depend on dDelta and so
    // never needs repeating), a reduced basis for expressing hybridization
    // changes as a combination of the Nbath_fine probe-energy poles
    // 1/(iw - e_[j]), then projects dDelta onto it into alphan_.
    // Equivalent to the original expand(complex*).
    void expand(const std::vector<Complex>& dDelta, Rng& rng);

    const Parameters* params_;
    const MatsubaraGrid* grid_;

    std::vector<double> t2u_raw_, t2d_raw_, eu_raw_, ed_raw_;
    double mu_raw_ = 0.0, h_raw_ = 0.0;
    std::vector<double> e_;  // Nbath_fine probe energies, fixed at construction

    std::vector<Complex> gu0_, gd0_;
    std::vector<Complex> gu_, gd_, gx_, ga_;

    std::vector<Complex> dg0_mu_, dg0_hz_, dgz_mu_, dgz_hz_, dgx_hx_, dga_ha_;
    // Each Nbath_fine x Nw: row j is the derivative w.r.t. the
    // hybridization pole at e_[j].
    Matrix<Complex> dg0_D0_, dg0_Dz_, dgz_D0_, dgz_Dz_, dgx_Dx_, dga_Da_;

    std::vector<Complex> alphan_;
    bool basis_ready_ = false;
    Matrix<double> pseudo1_;  // Nbath_fine x Nbath_fine

    // Reused scratch for get_g_static()/get_g(): two different Fock space
    // sizes (2*Nbath_raw+2 sites vs. 2*Nbath_raw+3, the extra site being
    // the probe bath level), each diagonalized dozens of times per
    // get_derivatives() call. As with RawImpuritySolver's
    // adjust_bath_couplings(), reusing one buffer per shape (rather than
    // allocating fresh per call) is a real performance requirement here,
    // not just style.
    int n_fock_static_, n_fock_probe_;
    Matrix<double> H_static_, psi_static_, cu_static_, cd_static_, cu_psi_static_, cd_psi_static_;
    std::vector<double> e_static_;
    Matrix<double> H_probe_, psi_probe_, cu_probe_, cd_probe_, cu_psi_probe_, cd_psi_probe_;
    std::vector<double> e_probe_;
};

}  // namespace dmft
