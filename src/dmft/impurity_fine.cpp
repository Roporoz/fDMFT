#include "impurity_fine.hpp"

#include <cmath>
#include <complex>

namespace dmft {

FineImpuritySolver::FineImpuritySolver(const Parameters& params, const MatsubaraGrid& grid)
    : g0(grid.size()),
      gz(grid.size()),
      gu(grid.size()),
      gd(grid.size()),
      gx(grid.size()),
      gy(grid.size()),
      gax(grid.size()),
      gay(grid.size()),
      Delta_u(grid.size(), 0.0),
      Delta_d(grid.size(), 0.0),
      Delta_x(grid.size(), 0.0),
      Delta_y(grid.size(), 0.0),
      Delta_ax(grid.size(), 0.0),
      Delta_ay(grid.size(), 0.0),
      params_(&params),
      grid_(&grid),
      t2u_raw_(params.Nbath_raw),
      t2d_raw_(params.Nbath_raw),
      eu_raw_(params.Nbath_raw),
      ed_raw_(params.Nbath_raw),
      e_(params.Nbath_fine),
      gu0_(grid.size()),
      gd0_(grid.size()),
      gu_(grid.size()),
      gd_(grid.size()),
      gx_(grid.size()),
      ga_(grid.size()),
      dg0_mu_(grid.size()),
      dg0_hz_(grid.size()),
      dgz_mu_(grid.size()),
      dgz_hz_(grid.size()),
      dgx_hx_(grid.size()),
      dga_ha_(grid.size()),
      dg0_D0_(params.Nbath_fine, grid.size()),
      dg0_Dz_(params.Nbath_fine, grid.size()),
      dgz_D0_(params.Nbath_fine, grid.size()),
      dgz_Dz_(params.Nbath_fine, grid.size()),
      dgx_Dx_(params.Nbath_fine, grid.size()),
      dga_Da_(params.Nbath_fine, grid.size()),
      alphan_(params.Nbath_fine),
      pseudo1_(params.Nbath_fine, params.Nbath_fine),
      n_fock_static_(1 << (2 * params.Nbath_raw + 2)),
      n_fock_probe_(1 << (2 * params.Nbath_raw + 3)),
      H_static_(n_fock_static_, n_fock_static_),
      psi_static_(n_fock_static_, n_fock_static_),
      cu_static_(n_fock_static_, n_fock_static_),
      cd_static_(n_fock_static_, n_fock_static_),
      cu_psi_static_(n_fock_static_, n_fock_static_),
      cd_psi_static_(n_fock_static_, n_fock_static_),
      e_static_(n_fock_static_, 0.0),
      H_probe_(n_fock_probe_, n_fock_probe_),
      psi_probe_(n_fock_probe_, n_fock_probe_),
      cu_probe_(n_fock_probe_, n_fock_probe_),
      cd_probe_(n_fock_probe_, n_fock_probe_),
      cu_psi_probe_(n_fock_probe_, n_fock_probe_),
      cd_psi_probe_(n_fock_probe_, n_fock_probe_),
      e_probe_(n_fock_probe_, 0.0) {
    for (int ne = 0; ne < params.Nbath_fine; ++ne)
        e_[ne] = params.Bath_fine_halfwidth * (-1.0 + ne * 2.0 / (params.Nbath_fine - 1.0));
}

void FineImpuritySolver::init(const RawImpuritySolver& raw, bool calc_g) {
    const int nb = params_->Nbath_raw;
    for (int l = 0; l < nb; ++l) {
        t2u_raw_[l] = raw.t2u[l];
        t2d_raw_[l] = raw.t2d[l];
        eu_raw_[l] = raw.eu[l];
        ed_raw_[l] = raw.ed[l];
    }
    mu_raw_ = raw.mu_loc;
    h_raw_ = raw.h_loc;

    const int n_freq = grid_->size();
    for (int w = 0; w < n_freq; ++w) {
        g0[w] = raw.gu[w] + raw.gd[w];
        gz[w] = raw.gu[w] - raw.gd[w];
        gu[w] = raw.gu[w];
        gd[w] = raw.gd[w];
        gx[w] = 0.0;
        gy[w] = 0.0;
        gax[w] = 0.0;
        gay[w] = 0.0;
    }

    if (calc_g) {
        get_g_static(0.0, false);
        gu0_ = gu_;
        gd0_ = gd_;
        get_derivatives();
    }
}

double FineImpuritySolver::magnetization() const {
    return extrapolated_matsubara_sum(grid_->size(),
                                       [&](int w) { return 2.0 * gz[w].real() / params_->beta; });
}

double FineImpuritySolver::occupation() const {
    return extrapolated_matsubara_sum(grid_->size(),
                                       [&](int w) { return 2.0 * g0[w].real() / params_->beta; });
}

void FineImpuritySolver::get_g_static(double dh, bool anomal) {
    const int nb = params_->Nbath_raw;
    const double U = params_->U;
    const double beta = params_->beta;
    const int n_fock = n_fock_static_;
    const int n_freq = grid_->size();

    for (int i = 0; i < n_fock; ++i)
        for (int j = 0; j < n_fock; ++j) {
            H_static_[i][j] = 0.0;
            cu_psi_static_[i][j] = 0.0;
            cd_psi_static_[i][j] = 0.0;
            cu_static_[i][j] = 0.0;
            cd_static_[i][j] = 0.0;
        }

    for (int j = 0; j < n_fock; ++j) {
        H_static_[j][j] += -mu_raw_ * (occ(j, 0) + occ(j, nb + 1)) -
                            h_raw_ * (occ(j, 0) - occ(j, nb + 1)) +
                            U * ((occ(j, 0) - 0.5) * (occ(j, nb + 1) - 0.5));

        if (anomal) {
            H_static_[j][j] += -dh * (occ(j, 0) - occ(j, nb + 1));
            if (occ(j, 0) == 0 && occ(j, nb + 1) == 0) {
                const int j1 = j + 1 + (1 << (nb + 1));
                int nu = 0;
                for (int l = 0; l < nb; ++l) nu += occ(j, l + 1);
                H_static_[j][j1] -= (1 - 2 * (nu % 2)) * dh;
                H_static_[j1][j] = H_static_[j][j1];
            }
        } else {
            H_static_[j][j] += -dh * (occ(j, 0) + occ(j, nb + 1));
            if (occ(j, 0) == 0 && occ(j, nb + 1) == 1) {
                const int j1 = j + 1 - (1 << (nb + 1));
                int nu = 0;
                for (int l = 0; l < nb; ++l) nu += occ(j, l + 1);
                H_static_[j][j1] -= (1 - 2 * (nu % 2)) * dh;
                H_static_[j1][j] = H_static_[j][j1];
            }
        }

        for (int l = 0; l < nb; ++l)
            H_static_[j][j] += occ(j, l + 1) * eu_raw_[l] + occ(j, nb + l + 2) * ed_raw_[l];

        for (int l = 0; l < nb; ++l) {
            int nu = 0, nd = 0;
            for (int l2 = 0; l2 < l; ++l2) {
                nu += occ(j, l2 + 1);
                nd += occ(j, nb + l2 + 2);
            }
            if (occ(j, 0) == 0 && occ(j, l + 1) == 1) {
                const int j1 = j + 1 - (1 << (l + 1));
                H_static_[j][j1] += std::sqrt(t2u_raw_[l]) * (1 - 2 * (nu % 2));
                H_static_[j1][j] = H_static_[j][j1];
            }
            if (occ(j, nb + 1) == 0 && occ(j, l + nb + 2) == 1) {
                const int j1 = j + (1 << (nb + 1)) - (1 << (l + nb + 2));
                H_static_[j][j1] += std::sqrt(t2d_raw_[l]) * (1 - 2 * (nd % 2));
                H_static_[j1][j] = H_static_[j][j1];
            }
        }
    }

    jacobi_eigensolver(H_static_, psi_static_, e_static_);

    double z = 0.0;
    for (int j = 0; j < n_fock; ++j) z += expl(-beta * e_static_[j]);

    for (int j = 0; j < n_fock; ++j) {
        if (occ(j, 0) == 1)
            for (int j1 = 0; j1 < n_fock; ++j1) cu_psi_static_[j1][j - 1] = psi_static_[j1][j];
        if (occ(j, nb + 1) == 1) {
            double f = 1.0;
            for (int l = 0; l <= nb; ++l) f *= 1 - 2 * occ(j, l);
            for (int j1 = 0; j1 < n_fock; ++j1)
                cd_psi_static_[j1][j - (1 << (nb + 1))] = f * psi_static_[j1][j];
        }
    }

    for (int w = 0; w < n_freq; ++w) {
        gu_[w] = 0.0;
        gd_[w] = 0.0;
        gx_[w] = 0.0;
        ga_[w] = 0.0;
    }

    for (int j1 = 0; j1 < n_fock; ++j1)
        for (int j2 = 0; j2 < n_fock; ++j2) {
            double d = 0.0;
            for (int j3 = 0; j3 < n_fock; ++j3) d += psi_static_[j1][j3] * cu_psi_static_[j2][j3];
            cu_static_[j1][j2] = d;
            d = 0.0;
            for (int j3 = 0; j3 < n_fock; ++j3) d += psi_static_[j1][j3] * cd_psi_static_[j2][j3];
            cd_static_[j1][j2] = d;
        }

    for (int j1 = 0; j1 < n_fock; ++j1)
        for (int j2 = 0; j2 < n_fock; ++j2) {
            const double w0 = (expl(-beta * e_static_[j1]) + expl(-beta * e_static_[j2])) / z;
            double x;
            x = cu_static_[j1][j2] * cu_static_[j1][j2] * w0;
            for (int w = 0; w < n_freq; ++w)
                gu_[w] += x / (grid_->iw(w) + e_static_[j1] - e_static_[j2]);

            x = cd_static_[j1][j2] * cd_static_[j1][j2] * w0;
            for (int w = 0; w < n_freq; ++w)
                gd_[w] += x / (grid_->iw(w) + e_static_[j1] - e_static_[j2]);

            x = cu_static_[j1][j2] * cd_static_[j1][j2] * w0;
            for (int w = 0; w < n_freq; ++w)
                gx_[w] += 2.0 * x / (grid_->iw(w) + e_static_[j1] - e_static_[j2]);

            x = cu_static_[j1][j2] * cd_static_[j2][j1] * w0;
            for (int w = 0; w < n_freq; ++w)
                ga_[w] += 2.0 * x / (grid_->iw(w) + e_static_[j1] - e_static_[j2]);
        }
}

void FineImpuritySolver::get_g(double t_, double e_probe, int flag) {
    const int nb = params_->Nbath_raw;
    const double U = params_->U;
    const double beta = params_->beta;
    const int n_fock = n_fock_probe_;
    const int n_freq = grid_->size();

    for (int i = 0; i < n_fock; ++i)
        for (int j = 0; j < n_fock; ++j) {
            H_probe_[i][j] = 0.0;
            cu_psi_probe_[i][j] = 0.0;
            cd_psi_probe_[i][j] = 0.0;
            cu_probe_[i][j] = 0.0;
            cd_probe_[i][j] = 0.0;
        }

    for (int j = 0; j < n_fock; ++j) {
        H_probe_[j][j] += -mu_raw_ * (occ(j, 0) + occ(j, nb + 1)) -
                           h_raw_ * (occ(j, 0) - occ(j, nb + 1)) +
                           U * ((occ(j, 0) - 0.5) * (occ(j, nb + 1) - 0.5));

        for (int l = 0; l < nb; ++l)
            H_probe_[j][j] += occ(j, l + 1) * eu_raw_[l] + occ(j, nb + l + 2) * ed_raw_[l];

        for (int l = 0; l < nb; ++l) {
            int nu = 0, nd = 0;
            for (int l2 = 0; l2 < l; ++l2) {
                nu += occ(j, l2 + 1);
                nd += occ(j, nb + l2 + 2);
            }
            if (occ(j, 0) == 0 && occ(j, l + 1) == 1) {
                const int j1 = j + 1 - (1 << (l + 1));
                H_probe_[j][j1] += std::sqrt(t2u_raw_[l]) * (1 - 2 * (nu % 2));
                H_probe_[j1][j] = H_probe_[j][j1];
            }
            if (occ(j, nb + 1) == 0 && occ(j, l + nb + 2) == 1) {
                const int j1 = j + (1 << (nb + 1)) - (1 << (l + nb + 2));
                H_probe_[j][j1] += std::sqrt(t2d_raw_[l]) * (1 - 2 * (nd % 2));
                H_probe_[j1][j] = H_probe_[j][j1];
            }
        }

        H_probe_[j][j] += e_probe * occ(j, 2 * nb + 2);
        if (flag == 0) {
            if (occ(j, 0) == 0 && occ(j, 2 * nb + 2) == 1) {
                const int j1 = j + 1 - (1 << (2 * nb + 2));
                int nu = 0;
                for (int l = 0; l < 2 * nb + 1; ++l) nu += occ(j, l + 1);
                H_probe_[j][j1] += (1 - 2 * (nu % 2)) * t_;
                H_probe_[j1][j] = H_probe_[j][j1];
            }
        }
        if (flag == 1) {
            if (occ(j, nb + 1) == 0 && occ(j, 2 * nb + 2) == 1) {
                const int j1 = j + (1 << (nb + 1)) - (1 << (2 * nb + 2));
                int nu = 0;
                for (int l = 0; l < nb; ++l) nu += occ(j, l + nb + 2);
                H_probe_[j][j1] += (1 - 2 * (nu % 2)) * t_;
                H_probe_[j1][j] = H_probe_[j][j1];
            }
        }
        if (flag == 2) {
            if (occ(j, 0) == 0 && occ(j, 2 * nb + 2) == 1) {
                const int j1 = j + 1 - (1 << (2 * nb + 2));
                int nu = 0;
                for (int l = 0; l < 2 * nb + 1; ++l) nu += occ(j, l + 1);
                H_probe_[j][j1] += (1 - 2 * (nu % 2)) * t_;
                H_probe_[j1][j] = H_probe_[j][j1];
            }
            if (occ(j, nb + 1) == 0 && occ(j, 2 * nb + 2) == 1) {
                const int j1 = j + (1 << (nb + 1)) - (1 << (2 * nb + 2));
                int nu = 0;
                for (int l = 0; l < nb; ++l) nu += occ(j, l + nb + 2);
                H_probe_[j][j1] += (1 - 2 * (nu % 2)) * t_;
                H_probe_[j1][j] = H_probe_[j][j1];
            }
            if (occ(j, nb + 1) == 1 && occ(j, 2 * nb + 2) == 1) {
                const int j1 = j - (1 << (nb + 1)) - (1 << (2 * nb + 2));
                int nu = 0;
                for (int l = 0; l < nb; ++l) nu += occ(j, l + nb + 2);
                H_probe_[j][j1] -= (1 - 2 * (nu % 2)) * t_;
                H_probe_[j1][j] = H_probe_[j][j1];
            }
        }
    }

    jacobi_eigensolver(H_probe_, psi_probe_, e_probe_);

    double z = 0.0;
    for (int j = 0; j < n_fock; ++j) z += expl(-beta * e_probe_[j]);

    for (int j = 0; j < n_fock; ++j) {
        if (occ(j, 0) == 1)
            for (int j1 = 0; j1 < n_fock; ++j1) cu_psi_probe_[j1][j - 1] = psi_probe_[j1][j];
        if (occ(j, nb + 1) == 1) {
            double f = 1.0;
            for (int l = 0; l <= nb; ++l) f *= 1 - 2 * occ(j, l);
            for (int j1 = 0; j1 < n_fock; ++j1)
                cd_psi_probe_[j1][j - (1 << (nb + 1))] = f * psi_probe_[j1][j];
        }
    }

    for (int w = 0; w < n_freq; ++w) {
        gu_[w] = 0.0;
        gd_[w] = 0.0;
        gx_[w] = 0.0;
        ga_[w] = 0.0;
    }

    for (int j1 = 0; j1 < n_fock; ++j1)
        for (int j2 = 0; j2 < n_fock; ++j2) {
            double d = 0.0;
            for (int j3 = 0; j3 < n_fock; ++j3) d += psi_probe_[j1][j3] * cu_psi_probe_[j2][j3];
            cu_probe_[j1][j2] = d;
            d = 0.0;
            for (int j3 = 0; j3 < n_fock; ++j3) d += psi_probe_[j1][j3] * cd_psi_probe_[j2][j3];
            cd_probe_[j1][j2] = d;
        }

    for (int j1 = 0; j1 < n_fock; ++j1)
        for (int j2 = 0; j2 < n_fock; ++j2) {
            const double w0 = (expl(-beta * e_probe_[j1]) + expl(-beta * e_probe_[j2])) / z;
            double x;
            x = cu_probe_[j1][j2] * cu_probe_[j1][j2] * w0;
            for (int w = 0; w < n_freq; ++w)
                gu_[w] += x / (grid_->iw(w) + e_probe_[j1] - e_probe_[j2]);

            x = cd_probe_[j1][j2] * cd_probe_[j1][j2] * w0;
            for (int w = 0; w < n_freq; ++w)
                gd_[w] += x / (grid_->iw(w) + e_probe_[j1] - e_probe_[j2]);

            x = cu_probe_[j1][j2] * cd_probe_[j1][j2] * w0;
            for (int w = 0; w < n_freq; ++w)
                gx_[w] += 2.0 * x / (grid_->iw(w) + e_probe_[j1] - e_probe_[j2]);

            x = cu_probe_[j1][j2] * cd_probe_[j2][j1] * w0;
            for (int w = 0; w < n_freq; ++w)
                ga_[w] += 2.0 * x / (grid_->iw(w) + e_probe_[j1] - e_probe_[j2]);
        }
}

void FineImpuritySolver::get_derivatives() {
    const int n_freq = grid_->size();
    const double dh = 0.001;

    get_g_static(0.5 * dh, false);
    for (int w = 0; w < n_freq; ++w) {
        dg0_mu_[w] = (gu_[w] + gd_[w]) / dh;
        dgz_mu_[w] = (gu_[w] - gd_[w]) / dh;
        dgx_hx_[w] = gx_[w] / dh;
    }
    get_g_static(-0.5 * dh, false);
    for (int w = 0; w < n_freq; ++w) {
        dg0_mu_[w] -= (gu_[w] + gd_[w]) / dh;
        dgz_mu_[w] -= (gu_[w] - gd_[w]) / dh;
        dgx_hx_[w] -= gx_[w] / dh;
    }

    get_g_static(0.5 * dh, true);
    for (int w = 0; w < n_freq; ++w) {
        dg0_hz_[w] = (gu_[w] + gd_[w]) / dh;
        dgz_hz_[w] = (gu_[w] - gd_[w]) / dh;
        dga_ha_[w] = ga_[w] / dh;
    }
    get_g_static(-0.5 * dh, true);
    for (int w = 0; w < n_freq; ++w) {
        dg0_hz_[w] -= (gu_[w] + gd_[w]) / dh;
        dgz_hz_[w] -= (gu_[w] - gd_[w]) / dh;
        dga_ha_[w] -= ga_[w] / dh;
    }

    const double t_ = 0.01;
    for (int ne = 0; ne < params_->Nbath_fine; ++ne) {
        get_g(t_, e_[ne], 0);
        for (int w = 0; w < n_freq; ++w) {
            dg0_D0_[ne][w] = 2.0 * (gu_[w] + gd_[w] - gu0_[w] - gd0_[w]) / (-2.0 * t_ * t_);
            dgz_D0_[ne][w] = 2.0 * (gu_[w] - gd_[w] - gu0_[w] + gd0_[w]) / (-2.0 * t_ * t_);
            dg0_Dz_[ne][w] = 2.0 * (gu_[w] + gd_[w] - gu0_[w] - gd0_[w]) / (-2.0 * t_ * t_);
            dgz_Dz_[ne][w] = 2.0 * (gu_[w] - gd_[w] - gu0_[w] + gd0_[w]) / (-2.0 * t_ * t_);
        }
        get_g(t_ * std::sqrt(2.0), e_[ne], 0);
        for (int w = 0; w < n_freq; ++w) {
            dg0_D0_[ne][w] -= 0.5 * (gu_[w] + gd_[w] - gu0_[w] - gd0_[w]) / (-2.0 * t_ * t_);
            dgz_D0_[ne][w] -= 0.5 * (gu_[w] - gd_[w] - gu0_[w] + gd0_[w]) / (-2.0 * t_ * t_);
            dg0_Dz_[ne][w] -= 0.5 * (gu_[w] + gd_[w] - gu0_[w] - gd0_[w]) / (-2.0 * t_ * t_);
            dgz_Dz_[ne][w] -= 0.5 * (gu_[w] - gd_[w] - gu0_[w] + gd0_[w]) / (-2.0 * t_ * t_);
        }

        get_g(t_, e_[ne], 1);
        for (int w = 0; w < n_freq; ++w) {
            dg0_D0_[ne][w] += 2.0 * (gu_[w] + gd_[w] - gu0_[w] - gd0_[w]) / (-2.0 * t_ * t_);
            dgz_D0_[ne][w] += 2.0 * (gu_[w] - gd_[w] - gu0_[w] + gd0_[w]) / (-2.0 * t_ * t_);
            dg0_Dz_[ne][w] -= 2.0 * (gu_[w] + gd_[w] - gu0_[w] - gd0_[w]) / (-2.0 * t_ * t_);
            dgz_Dz_[ne][w] -= 2.0 * (gu_[w] - gd_[w] - gu0_[w] + gd0_[w]) / (-2.0 * t_ * t_);
        }
        get_g(t_ * std::sqrt(2.0), e_[ne], 1);
        for (int w = 0; w < n_freq; ++w) {
            dg0_D0_[ne][w] -= 0.5 * (gu_[w] + gd_[w] - gu0_[w] - gd0_[w]) / (-2.0 * t_ * t_);
            dgz_D0_[ne][w] -= 0.5 * (gu_[w] - gd_[w] - gu0_[w] + gd0_[w]) / (-2.0 * t_ * t_);
            dg0_Dz_[ne][w] += 0.5 * (gu_[w] + gd_[w] - gu0_[w] - gd0_[w]) / (-2.0 * t_ * t_);
            dgz_Dz_[ne][w] += 0.5 * (gu_[w] - gd_[w] - gu0_[w] + gd0_[w]) / (-2.0 * t_ * t_);
        }

        get_g(t_, e_[ne], 2);
        for (int w = 0; w < n_freq; ++w) {
            dgx_Dx_[ne][w] = 2.0 * gx_[w] / (-2.0 * t_ * t_);
            dga_Da_[ne][w] = 2.0 * ga_[w] / (-2.0 * t_ * t_);
        }
        get_g(t_ * std::sqrt(2.0), e_[ne], 2);
        for (int w = 0; w < n_freq; ++w) {
            dgx_Dx_[ne][w] -= 0.5 * gx_[w] / (-2.0 * t_ * t_);
            dga_Da_[ne][w] -= 0.5 * ga_[w] / (-2.0 * t_ * t_);
        }
    }
}

double FineImpuritySolver::scalar(const std::vector<Complex>& r1, const std::vector<Complex>& r2) const {
    return 2.0 * extrapolated_matsubara_sum(grid_->size(), [&](int w) {
                      return std::conj(r1[w]) * r2[w];
                  }).real() /
           params_->beta;
}

void FineImpuritySolver::expand(const std::vector<Complex>& dDelta, Rng& rng) {
    const int nf = params_->Nbath_fine;
    const int n_freq = grid_->size();

    if (!basis_ready_) {
        Matrix<double> m(nf, nf);
        std::vector<Complex> r1(n_freq), r2(n_freq);
        for (int i = 0; i < nf; ++i)
            for (int j = 0; j < nf; ++j) {
                for (int w = 0; w < n_freq; ++w) {
                    r1[w] = 1.0 / (grid_->iw(w) - e_[i]);
                    r2[w] = 1.0 / (grid_->iw(w) - e_[j]);
                }
                m[i][j] = scalar(r1, r2);
            }

        Matrix<double> psik(nf, nf);
        std::vector<double> e_eff(nf);
        const int n_k = leading_eigenvectors(m, psik, e_eff, nf, rng, 1e-5);

        for (int k = 0; k < n_k; ++k) {
            for (int w = 0; w < n_freq; ++w) {
                r1[w] = 0.0;
                for (int j = 0; j < nf; ++j) r1[w] += psik[k][j] / (grid_->iw(w) - e_[j]);
            }
            e_eff[k] = scalar(r1, r1);
        }

        for (int i = 0; i < nf; ++i)
            for (int j = 0; j < nf; ++j) {
                double p = 0.0;
                for (int k = 0; k < n_k; ++k) p += psik[k][i] * psik[k][j] / e_eff[k];
                pseudo1_[i][j] = p;
            }

        basis_ready_ = true;
    }

    for (int j = 0; j < nf; ++j) alphan_[j] = 0.0;

    std::vector<Complex> r1(n_freq);
    for (int j_ = 0; j_ < nf; ++j_) {
        for (int w = 0; w < n_freq; ++w) r1[w] = -1.0 / (grid_->iw(w) - e_[j_]);
        const double r = scalar(r1, dDelta);
        for (int j = 0; j < nf; ++j) alphan_[j] += pseudo1_[j][j_] * r;
    }
}

void FineImpuritySolver::set_g_raw() {
    const int n_freq = grid_->size();
    for (int w = 0; w < n_freq; ++w) {
        g0[w] = gu0_[w] + gd0_[w];
        gz[w] = gu0_[w] - gd0_[w];
        gu[w] = gu0_[w];
        gd[w] = gd0_[w];
        gx[w] = 0.0;
        gy[w] = 0.0;
        gax[w] = 0.0;
        gay[w] = 0.0;
    }
}

void FineImpuritySolver::set_g(Rng& rng) {
    const int nb = params_->Nbath_raw;
    const int nf = params_->Nbath_fine;
    const int n_freq = grid_->size();

    for (int w = 0; w < n_freq; ++w) {
        g0[w] = gu0_[w] + gd0_[w];
        gz[w] = gu0_[w] - gd0_[w];
        gu[w] = gu0_[w];
        gd[w] = gd0_[w];
        gx[w] = 0.0;
        gy[w] = 0.0;
        gax[w] = 0.0;
        gay[w] = 0.0;
    }

    std::vector<Complex> dDelta_0(n_freq), dDelta_z(n_freq);
    for (int w = 0; w < n_freq; ++w) {
        dDelta_0[w] = Delta_u[w] + Delta_d[w];
        dDelta_z[w] = Delta_u[w] - Delta_d[w];
        for (int j = 0; j < nb; ++j) {
            dDelta_0[w] -= t2u_raw_[j] / (grid_->iw(w) - eu_raw_[j]) +
                           t2d_raw_[j] / (grid_->iw(w) - ed_raw_[j]);
            dDelta_z[w] -= t2u_raw_[j] / (grid_->iw(w) - eu_raw_[j]) -
                           t2d_raw_[j] / (grid_->iw(w) - ed_raw_[j]);
        }
    }

    expand(dDelta_0, rng);
    for (int w = 0; w < n_freq; ++w) {
        for (int j = 0; j < nf; ++j) g0[w] += dg0_D0_[j][w] * alphan_[j];
        g0[w] += dg0_mu_[w] * (mu_fine - mu_raw_);
    }
    for (int w = 0; w < n_freq; ++w) {
        for (int j = 0; j < nf; ++j) gz[w] += dgz_D0_[j][w] * alphan_[j];
        gz[w] += dgz_mu_[w] * (mu_fine - mu_raw_);
    }

    expand(dDelta_z, rng);
    for (int w = 0; w < n_freq; ++w) {
        for (int j = 0; j < nf; ++j) g0[w] += dg0_Dz_[j][w] * alphan_[j];
        g0[w] += dg0_hz_[w] * (hz_fine - h_raw_);
    }
    for (int w = 0; w < n_freq; ++w) {
        for (int j = 0; j < nf; ++j) gz[w] += dgz_Dz_[j][w] * alphan_[j];
        gz[w] += dgz_hz_[w] * (hz_fine - h_raw_);
    }

    for (int w = 0; w < n_freq; ++w) {
        gu[w] = 0.5 * (g0[w] + gz[w]);
        gd[w] = 0.5 * (g0[w] - gz[w]);
    }

    expand(Delta_x, rng);
    for (int w = 0; w < n_freq; ++w) {
        for (int j = 0; j < nf; ++j) gx[w] += dgx_Dx_[j][w] * alphan_[j];
        gx[w] += dgx_hx_[w] * hx_fine;
    }
    expand(Delta_y, rng);
    for (int w = 0; w < n_freq; ++w) {
        for (int j = 0; j < nf; ++j) gy[w] += dgx_Dx_[j][w] * alphan_[j];
        gy[w] += dgx_hx_[w] * hx_fine;
    }

    expand(Delta_ax, rng);
    for (int w = 0; w < n_freq; ++w) {
        for (int j = 0; j < nf; ++j) gax[w] += dga_Da_[j][w] * alphan_[j];
        gax[w] += dga_ha_[w] * hax_fine;
    }
    expand(Delta_ay, rng);
    for (int w = 0; w < n_freq; ++w) {
        for (int j = 0; j < nf; ++j) gay[w] += dga_Da_[j][w] * alphan_[j];
        gay[w] += dga_ha_[w] * hay_fine;
    }
}

}  // namespace dmft
