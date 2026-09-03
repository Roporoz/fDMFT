#include "impurity_raw.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <utility>

namespace dmft {

RawImpuritySolver::RawImpuritySolver(const Parameters& params, const MatsubaraGrid& grid)
    : t2u(params.Nbath_raw, 0.0),
      t2d(params.Nbath_raw, 0.0),
      eu(params.Nbath_raw),
      ed(params.Nbath_raw),
      gu(grid.size()),
      gd(grid.size()),
      sigmau(grid.size(), 0.0),
      sigmad(grid.size(), 0.0),
      gudh(grid.size(), 0.0),
      gddh(grid.size(), 0.0),
      gxdh(grid.size(), 0.0),
      params_(&params),
      grid_(&grid),
      n_fock_(1 << (2 * params.Nbath_raw + 2)),
      H_(n_fock_, n_fock_),
      Psi_(n_fock_, n_fock_),
      cu_psi_(n_fock_, n_fock_),
      cd_psi_(n_fock_, n_fock_),
      cu_(n_fock_, n_fock_),
      cd_(n_fock_, n_fock_),
      E_(n_fock_, 0.0),
      nud_(n_fock_, 0),
      bath_fit_matrix_(params.Nbath_raw, params.Nbath_raw),
      newton_hessian_(params.Nbath_raw, params.Nbath_raw),
      newton_gradient_(params.Nbath_raw) {
    for (int w = 0; w < grid.size(); ++w) {
        gu[w] = 1.0 / grid.iw(w);
        gd[w] = 1.0 / grid.iw(w);
    }
    for (int l = 0; l < params.Nbath_raw; ++l) {
        eu[l] = 0.1 * l;
        ed[l] = 0.1 * l;
    }
}

double RawImpuritySolver::init(const std::vector<Complex>& delta_up,
                                const std::vector<Complex>& delta_down, double mu_loc_,
                                double h_loc_, Rng& rng) {
    tolerance = fit_bath(delta_up, delta_down, rng);
    return diagonalize(mu_loc_, h_loc_);
}

double RawImpuritySolver::fit_bath(const std::vector<Complex>& delta_up,
                                    const std::vector<Complex>& delta_down, Rng& rng) {
    return adjust_bath_levels(delta_up, t2u, eu, rng) + adjust_bath_levels(delta_down, t2d, ed, rng);
}

double RawImpuritySolver::adjust_bath_couplings(const std::vector<Complex>& delta_w,
                                                 std::vector<double>& t2_bath,
                                                 const std::vector<double>& e_bath) const {
    const int nb = params_->Nbath_raw;
    const int n_ed = params_->NwED;

    Matrix<double>& A = bath_fit_matrix_;
    for (int i = 0; i < nb; ++i)
        for (int j = i; j < nb; ++j) {
            const double a = extrapolated_matsubara_sum(n_ed, [&](int w) {
                return (e_bath[i] / (sqr(e_bath[i]) - grid_->iw2(w)) +
                        e_bath[j] / (sqr(e_bath[j]) - grid_->iw2(w))) /
                       (e_bath[i] + e_bath[j]);
            });
            A[i][j] = a;
            A[j][i] = a;
        }

    t2_bath.resize(nb);
    for (int i = 0; i < nb; ++i) {
        t2_bath[i] = extrapolated_matsubara_sum(n_ed, [&](int w) {
            return std::real(std::conj(delta_w[w]) / (grid_->iw(w) - e_bath[i]));
        });
    }

    if (!solve_in_place(A, t2_bath)) return 100.0;

    const double chi2 = extrapolated_matsubara_sum(n_ed, [&](int w) {
        Complex r = delta_w[w];
        for (int i = 0; i < nb; ++i) r -= t2_bath[i] / (grid_->iw(w) - e_bath[i]);
        return norm2(r);
    });
    return 2.0 * chi2 / params_->beta;
}

double RawImpuritySolver::adjust_bath_levels(const std::vector<Complex>& delta_w,
                                              std::vector<double>& t2_bath,
                                              std::vector<double>& e_bath, Rng& rng) {
    const int nb = params_->Nbath_raw;
    const std::vector<double> e_prev = e_bath;

    auto randomize = [&] {
        if (rng.uniform01() < 0.5) {
            for (int j = 0; j < nb; ++j) e_bath[j] = jump_scale_ * rng.gaussian2d().real();
        } else {
            for (int j = 0; j < nb; ++j)
                e_bath[j] = e_prev[j] + 0.3 * jump_scale_ * rng.gaussian2d().real();
        }
    };

    randomize();
    double s_min = adjust_bath_couplings(delta_w, t2_bath, e_bath);
    std::vector<double> e_bath_min = e_bath;

    int i_min = -1;
    bool t2_positive = false;
    double s = s_min;

    for (int i = 0; i < trial_budget_; ++i) {
        randomize();
        s = adjust_bath_couplings(delta_w, t2_bath, e_bath);
        if (!(s < 0) && s < 0.9999 * s_min) {
            s_min = s;
            e_bath_min = e_bath;
            i_min = i;
            t2_positive = true;
            for (int l = 0; l < nb; ++l) t2_positive = t2_positive && t2_bath[l] > 0;
        }

        // Newton refinement: locally minimize the chi^2 residual s(e_bath)
        // via a finite-difference gradient/Hessian, for up to 50 steps or
        // until it stops improving.
        const double fd_step = 0.001;
        std::vector<double>& grad = newton_gradient_;
        Matrix<double>& hess = newton_hessian_;
        for (int ii = 0; ii < 50; ++ii) {
            for (int j = 0; j < nb; ++j) {
                e_bath[j] += fd_step;
                const double sp = adjust_bath_couplings(delta_w, t2_bath, e_bath);
                e_bath[j] -= 2.0 * fd_step;
                const double sm = adjust_bath_couplings(delta_w, t2_bath, e_bath);
                e_bath[j] += fd_step;
                grad[j] = (sp - sm) / (2.0 * fd_step);
                hess[j][j] = (sp + sm - 2.0 * s) / (fd_step * fd_step);
            }
            for (int j1 = 0; j1 < nb; ++j1)
                for (int j2 = j1 + 1; j2 < nb; ++j2) {
                    e_bath[j1] += fd_step;
                    e_bath[j2] += fd_step;
                    const double sxy = adjust_bath_couplings(delta_w, t2_bath, e_bath);
                    e_bath[j1] -= fd_step;
                    e_bath[j2] -= fd_step;

                    e_bath[j1] -= fd_step;
                    e_bath[j2] -= fd_step;
                    const double sxym = adjust_bath_couplings(delta_w, t2_bath, e_bath);
                    e_bath[j1] += fd_step;
                    e_bath[j2] += fd_step;

                    hess[j1][j2] = (sxy + sxym - 2.0 * s - hess[j1][j1] * fd_step * fd_step -
                                    hess[j2][j2] * fd_step * fd_step) /
                                   (2.0 * fd_step * fd_step);
                    hess[j2][j1] = hess[j1][j2];
                }

            if (solve_in_place(hess, grad))
                for (int j = 0; j < nb; ++j) e_bath[j] -= grad[j];
            const double s0 = s;
            s = adjust_bath_couplings(delta_w, t2_bath, e_bath);
            if (s >= s0) break;

            if (!(s < 0) && s < 0.9999 * s_min) {
                s_min = s;
                e_bath_min = e_bath;
                i_min = i;
                t2_positive = true;
                for (int l = 0; l < nb; ++l) t2_positive = t2_positive && t2_bath[l] > 0;
            }
        }

        if (s_min < 1e-10 && t2_positive) break;
        if (!(s < 0) && s < 1e-10 && t2_positive) {
            s_min = s;
            e_bath_min = e_bath;
            i_min = i;
            break;
        }
    }

    e_bath = e_bath_min;
    const double chi2 = adjust_bath_couplings(delta_w, t2_bath, e_bath);
    if (!t2_positive) {
        std::cerr << "warning: adjust_bath_levels landed on an unphysical bath\n";
        for (int l = 0; l < nb; ++l) std::cerr << "  e=" << e_bath[l] << "  t2=" << t2_bath[l] << "\n";
    }

    double jump_max = 0.0;
    for (int j = 0; j < nb; ++j) jump_max = std::max(jump_max, std::abs(e_bath_min[j]));
    if (i_min >= 0) trial_budget_ = static_cast<int>(0.9 * trial_budget_ + i_min);
    jump_scale_ = 0.9 * jump_scale_ + 0.1 * jump_max;
    if (trial_budget_ < 300) trial_budget_ = 300;

    return chi2;
}

double RawImpuritySolver::diagonalize(double mu_loc_, double h_loc_) {
    mu_loc = mu_loc_;
    h_loc = h_loc_;

    const int nb = params_->Nbath_raw;
    const double U = params_->U;
    const double beta = params_->beta;

    for (int i = 0; i < n_fock_; ++i)
        for (int j = 0; j < n_fock_; ++j) {
            H_[i][j] = 0.0;
            cu_psi_[i][j] = 0.0;
            cd_psi_[i][j] = 0.0;
            cu_[i][j] = 0.0;
            cd_[i][j] = 0.0;
        }
    for (int j = 0; j < n_fock_; ++j) {
        nud_[j] = 0;
        for (int l = 0; l <= nb; ++l) nud_[j] += occ(j, l) + (nb + 2) * occ(j, l + nb + 1);
    }

    for (int j = 0; j < n_fock_; ++j) {
        H_[j][j] += -mu_loc * (occ(j, 0) + occ(j, nb + 1)) - h_loc * (occ(j, 0) - occ(j, nb + 1)) +
                    U * ((occ(j, 0) - 0.5) * (occ(j, nb + 1) - 0.5));
        for (int l = 0; l < nb; ++l) H_[j][j] += occ(j, l + 1) * eu[l] + occ(j, nb + l + 2) * ed[l];
        for (int l = 0; l < nb; ++l) {
            int nu = 0, nd = 0;
            for (int l2 = 0; l2 < l; ++l2) {
                nu += occ(j, l2 + 1);
                nd += occ(j, nb + l2 + 2);
            }
            if (occ(j, 0) == 0 && occ(j, l + 1) == 1) {
                const int j1 = j + 1 - (1 << (l + 1));
                H_[j][j1] += std::sqrt(t2u[l]) * (1 - 2 * (nu % 2));
                H_[j1][j] = H_[j][j1];
            }
            if (occ(j, nb + 1) == 0 && occ(j, l + nb + 2) == 1) {
                const int j1 = j + (1 << (nb + 1)) - (1 << (l + nb + 2));
                H_[j][j1] += std::sqrt(t2d[l]) * (1 - 2 * (nd % 2));
                H_[j1][j] = H_[j][j1];
            }
        }
    }

    Matrix<double> psi;
    std::vector<double> energies;
    block_eigensolver(H_, psi, energies, nud_, sqr(nb + 2));
    Psi_ = std::move(psi);
    E_ = std::move(energies);

    Z_ = 0.0;
    for (int j = 0; j < n_fock_; ++j) Z_ += expl(-beta * E_[j]);

    long double lnZ = logl(Z_);
    for (int l = 0; l < nb; ++l)
        lnZ -= std::log(1.0 + std::exp(-beta * eu[l])) + std::log(1.0 + std::exp(-beta * ed[l]));

    for (int j = 0; j < n_fock_; ++j) {
        if (occ(j, 0) == 1)
            for (int j1 = 0; j1 < n_fock_; ++j1) cu_psi_[j1][j - 1] = Psi_[j1][j];
        if (occ(j, nb + 1) == 1) {
            const double f = 1 - 2 * ((nud_[j] % (nb + 2)) % 2);
            for (int j1 = 0; j1 < n_fock_; ++j1) cd_psi_[j1][j - (1 << (nb + 1))] = f * Psi_[j1][j];
        }
    }

    const int n_freq = grid_->size();
    gu.assign(n_freq, 0.0);
    gd.assign(n_freq, 0.0);

    for (int j1 = 0; j1 < n_fock_; ++j1)
        for (int j2 = 0; j2 < n_fock_; ++j2) {
            if (nud_[j1] == nud_[j2] - 1) {
                double d = 0.0;
                for (int j3 = 0; j3 < n_fock_; ++j3) d += Psi_[j1][j3] * cu_psi_[j2][j3];
                cu_[j1][j2] = d;
                const double x = d * d * ((expl(-beta * E_[j1]) + expl(-beta * E_[j2])) / Z_);
                for (int w = 0; w < n_freq; ++w) gu[w] += x / (grid_->iw(w) + E_[j1] - E_[j2]);
            }
            if (nud_[j1] == nud_[j2] - (nb + 2)) {
                double d = 0.0;
                for (int j3 = 0; j3 < n_fock_; ++j3) d += Psi_[j1][j3] * cd_psi_[j2][j3];
                cd_[j1][j2] = d;
                const double x = d * d * ((expl(-beta * E_[j1]) + expl(-beta * E_[j2])) / Z_);
                for (int w = 0; w < n_freq; ++w) gd[w] += x / (grid_->iw(w) + E_[j1] - E_[j2]);
            }
        }

    sigmau.assign(n_freq, 0.0);
    sigmad.assign(n_freq, 0.0);
    for (int w = 0; w < n_freq; ++w) {
        sigmau[w] = grid_->iw(w) - 1.0 / gu[w] + h_loc + mu_loc;
        sigmad[w] = grid_->iw(w) - 1.0 / gd[w] - h_loc + mu_loc;
        for (int l = 0; l < nb; ++l) {
            sigmau[w] -= t2u[l] / (grid_->iw(w) - eu[l]);
            sigmad[w] -= t2d[l] / (grid_->iw(w) - ed[l]);
        }
    }

    lnZ -= beta * mu_loc;
    return static_cast<double>(lnZ);
}

double RawImpuritySolver::magnetization() const {
    return extrapolated_matsubara_sum(grid_->size(), [&](int w) {
        return 2.0 * (gu[w] - gd[w]).real() / params_->beta;
    });
}

double RawImpuritySolver::occupation() const {
    return extrapolated_matsubara_sum(grid_->size(), [&](int w) {
        return 2.0 * (gu[w] + gd[w]).real() / params_->beta;
    });
}

double RawImpuritySolver::dm_dh_z() const {
    return extrapolated_matsubara_sum(grid_->size(), [&](int w) {
        return 2.0 * (gudh[w] - gddh[w]).real() / params_->beta;
    });
}

double RawImpuritySolver::dm_dh_x() const {
    return extrapolated_matsubara_sum(
        grid_->size(), [&](int w) { return 2.0 * gxdh[w].real() / params_->beta; });
}

void RawImpuritySolver::store(std::vector<Complex>& delta_u_out, std::vector<Complex>& delta_d_out,
                               std::vector<Complex>& g_u_out, std::vector<Complex>& g_d_out) const {
    const int n_freq = grid_->size();
    delta_u_out.assign(n_freq, 0.0);
    delta_d_out.assign(n_freq, 0.0);
    for (int w = 0; w < n_freq; ++w) {
        for (int l = 0; l < params_->Nbath_raw; ++l) {
            delta_u_out[w] += t2u[l] / (grid_->iw(w) - eu[l]);
            delta_d_out[w] += t2d[l] / (grid_->iw(w) - ed[l]);
        }
    }
    g_u_out = gu;
    g_d_out = gd;
}

}  // namespace dmft
