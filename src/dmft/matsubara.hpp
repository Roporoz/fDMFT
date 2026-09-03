#pragma once

#include <vector>

#include "core.hpp"

// Replaces legacy solver.cpp's global Iw/Iw2 arrays and set_Iw(), which
// recomputed both from scratch on every call to solver_raw::init(),
// solver_fine::get_g_static() and solver_fine::get_g() -- i.e. on every
// single exact diagonalization, needlessly, since Nw and beta don't change
// between those calls. Here it's built once and passed around.
namespace dmft {

class MatsubaraGrid {
public:
    MatsubaraGrid(int n_freq, double beta) : iw_(n_freq), iw2_(n_freq) {
        for (int w = 0; w < n_freq; ++w) {
            const double wn = (2.0 * w + 1) * pi / beta;
            iw_[w] = I * wn;
            iw2_[w] = -wn * wn;
        }
    }

    int size() const { return static_cast<int>(iw_.size()); }
    Complex iw(int w) const { return iw_[w]; }
    double iw2(int w) const { return iw2_[w]; }

private:
    std::vector<Complex> iw_;
    std::vector<double> iw2_;
};

// The original repeats this exact pattern everywhere it needs a Matsubara
// sum extrapolated to Nw -> infinity: accumulate `term(w)` over
// w in [0, n_terms), remembering the running total at the Nw/2 and 3Nw/4
// marks, then combine all three as 8*a - 9*a3 + 2*a2 (a 3-point Richardson-
// style extrapolation). Appeared verbatim, independently, in solver_raw::
// s()/n()/dszdh()/dsxdh()/adjust_bath_couplings()/adjust_bath_levels() --
// this factors that one pattern out instead of re-deriving it each time.
template <typename Term>
double extrapolated_matsubara_sum(int n_terms, Term term) {
    double a = 0.0, a2 = 0.0, a3 = 0.0;
    for (int w = 0; w < n_terms; ++w) {
        a += term(w);
        if (w + 1 == n_terms / 2) a2 = a;
        if (w + 1 == 3 * n_terms / 4) a3 = a;
    }
    return 8.0 * a - 9.0 * a3 + 2.0 * a2;
}

}  // namespace dmft
