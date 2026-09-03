#include "dmft_loop.hpp"

#include <cmath>
#include <complex>

#include "lattice.hpp"

namespace dmft {

double bath_fine_to_raw(const FineImpuritySolver& fine, RawImpuritySolver& raw, Rng& rng,
                         double lambda0) {
    return raw.init(fine.Delta_u, fine.Delta_d, fine.mu_fine,
                     fine.hz_fine - lambda0 * fine.magnetization(), rng);
}

void bath_raw_to_fine(const RawImpuritySolver& raw, FineImpuritySolver& fine, Rng& rng) {
    fine.init(raw);
    fine.set_g(rng);
}

AntiferromagneticDMFTLoop::AntiferromagneticDMFTLoop(const Parameters& params,
                                                       const MatsubaraGrid& grid, double lambda0)
    : Sr(params, grid), Sf(params, grid), params_(&params), grid_(&grid), lambda0_(lambda0) {}

void AntiferromagneticDMFTLoop::init(double hz, Rng& rng) {
    for (int w = 0; w < grid_->size(); ++w) {
        // [sic] both lines assign Delta_u; Delta_d is never set here and
        // stays at the zero FineImpuritySolver's constructor gives it.
        // Preserved from the original bath_f2r-era solver.cpp verbatim
        // rather than silently fixed, since it's this loop's starting
        // guess and changing it changes the iteration's path (not
        // necessarily its converged fixed point, but that's a physics
        // question, not a porting one).
        Sf.Delta_u[w] = 0.1 / (grid_->iw(w) + 0.1 * I);
        Sf.Delta_u[w] = 0.1 / (grid_->iw(w) + 0.1 * I);
    }
    Sf.mu_fine = params_->mu;
    Sf.hz_fine = hz;
    bath_fine_to_raw(Sf, Sr, rng, lambda0_);
}

void AntiferromagneticDMFTLoop::iter_Delta(double factor) {
    const int n_freq = grid_->size();
    const int Lx = params_->Lx;
    const int Ly = params_->Ly;
    const double N = static_cast<double>(params_->N());

    for (int w = 0; w < n_freq; ++w) {
        Complex su = 0.0, sd = 0.0;
        for (int kx = 0; kx < Lx; ++kx)
            for (int ky = 0; ky < Ly / 2; ++ky) {
                const Complex au = 1.0 / Sf.gu[w] + Sf.Delta_u[w];
                const Complex ad = 1.0 / Sf.gd[w] + Sf.Delta_d[w];
                const Complex g00 = -dispersion(*params_, kx, ky) + 0.5 * (au + ad);
                const Complex g11 =
                    -dispersion(*params_, kx + Lx / 2, ky + Ly / 2) + 0.5 * (au + ad);
                const Complex g01 = 0.5 * (au - ad);
                const Complex g10 = g01;
                const Complex d = g00 * g11 - g01 * g10;
                su += (g11 + g00 - 2.0 * g10) / (d * N);
                sd += (g11 + g00 + 2.0 * g10) / (d * N);
            }
        Sf.Delta_u[w] += factor * (1.0 / Sf.gu[w] - 1.0 / su);
        Sf.Delta_d[w] += factor * (1.0 / Sf.gd[w] - 1.0 / sd);
        if (!af_flag) {
            Sf.Delta_u[w] = 0.5 * (Sf.Delta_u[w] + Sf.Delta_d[w]);
            Sf.Delta_d[w] = Sf.Delta_u[w];
        }
    }
}

int AntiferromagneticDMFTLoop::raw_loops(double factor, Rng& rng) {
    Complex g_prev = -1.0;
    for (int i = 0; i < 1000; ++i) {
        bath_fine_to_raw(Sf, Sr, rng, lambda0_);
        Sf.init(Sr, false);
        iter_Delta(factor);
        if (std::abs(Sf.gu[0] - g_prev) > 1e-4) {
            g_prev = Sf.gu[0];
        } else {
            return i;
        }
    }
    return -1;
}

int AntiferromagneticDMFTLoop::fine_loops(double factor, Rng& rng) {
    Complex g_prev = -1.0;
    Sf.init(Sr);
    const double h_prev = Sf.hz_fine;
    for (int i = 0; i < 1000; ++i) {
        Sf.set_g(rng);
        iter_Delta(factor);
        Sf.hz_fine = h_prev - lambda0_ * Sf.magnetization();
        if (std::abs(Sf.gu[0] - g_prev) > 1e-6) {
            g_prev = Sf.gu[0];
        } else {
            Sf.hz_fine = h_prev;
            return i;
        }
    }
    Sf.hz_fine = h_prev;
    return -1;
}

}  // namespace dmft
