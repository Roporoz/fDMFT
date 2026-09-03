// Cross-checks the src/dmft/linalg.hpp rewrite against the original
// headers.new.h implementations (via legacy_bridge.cpp) on identical
// inputs. See legacy/README or the commit that introduced this file for
// why this replaces a whole-binary before/after comparison: the original
// program seeds its RNG from time(NULL), so two runs of the unmodified
// binary already disagree with each other.
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

#include "../src/dmft/linalg.hpp"
#include "legacy_bridge.hpp"

namespace {

using dmft::Complex;
using dmft::Matrix;

int g_failures = 0;

void expect(bool condition, const char* what) {
    if (!condition) {
        std::printf("FAIL: %s\n", what);
        ++g_failures;
    } else {
        std::printf("ok:   %s\n", what);
    }
}

double max_abs_diff(const std::vector<double>& a, const std::vector<double>& b) {
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
    return m;
}

double max_abs_diff(const std::vector<Complex>& a, const std::vector<Complex>& b) {
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
    return m;
}

std::vector<double> flatten(const Matrix<double>& m) {
    std::vector<double> flat(static_cast<std::size_t>(m.n_rows()) * m.n_cols());
    for (int i = 0; i < m.n_rows(); ++i)
        for (int j = 0; j < m.n_cols(); ++j) flat[i * m.n_cols() + j] = m[i][j];
    return flat;
}

std::vector<Complex> flatten(const Matrix<Complex>& m) {
    std::vector<Complex> flat(static_cast<std::size_t>(m.n_rows()) * m.n_cols());
    for (int i = 0; i < m.n_rows(); ++i)
        for (int j = 0; j < m.n_cols(); ++j) flat[i * m.n_cols() + j] = m[i][j];
    return flat;
}

Matrix<double> random_matrix(dmft::Rng& rng, int n) {
    Matrix<double> m(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) m[i][j] = 2.0 * rng.uniform01() - 1.0;
    return m;
}

Matrix<Complex> random_complex_matrix(dmft::Rng& rng, int n) {
    Matrix<Complex> m(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            m[i][j] = Complex(2.0 * rng.uniform01() - 1.0, 2.0 * rng.uniform01() - 1.0);
    return m;
}

Matrix<double> random_symmetric_matrix(dmft::Rng& rng, int n) {
    Matrix<double> m(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j) {
            const double v = 2.0 * rng.uniform01() - 1.0;
            m[i][j] = v;
            m[j][i] = v;
        }
    return m;
}

std::vector<double> random_vector(dmft::Rng& rng, int n) {
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i) v[i] = 2.0 * rng.uniform01() - 1.0;
    return v;
}

std::vector<Complex> random_complex_vector(dmft::Rng& rng, int n) {
    std::vector<Complex> v(n);
    for (int i = 0; i < n; ++i) v[i] = Complex(2.0 * rng.uniform01() - 1.0, 2.0 * rng.uniform01() - 1.0);
    return v;
}

void test_invert() {
    const int n = 6;
    dmft::Rng rng(1);
    Matrix<Complex> a = random_complex_matrix(rng, n);

    Matrix<Complex> a_new = a;
    dmft::invert(a_new);

    std::vector<Complex> a_legacy_flat = flatten(a);
    legacy_bridge::invert(a_legacy_flat, n);

    expect(max_abs_diff(flatten(a_new), a_legacy_flat) < 1e-9, "invert matches legacy Inverse()");
}

void test_solve_complex() {
    const int n = 6;
    dmft::Rng rng(2);
    Matrix<Complex> a = random_complex_matrix(rng, n);
    std::vector<Complex> r = random_complex_vector(rng, n);

    Matrix<Complex> a_new = a;
    std::vector<Complex> x_new = r;
    dmft::solve_in_place(a_new, x_new);

    std::vector<Complex> a_legacy_flat = flatten(a);
    std::vector<Complex> x_legacy = r;
    legacy_bridge::solve_in_place_complex(a_legacy_flat, x_legacy, n);

    expect(max_abs_diff(x_new, x_legacy) < 1e-9, "solve_in_place<Complex> matches legacy div_left()");
}

void test_solve_double() {
    const int n = 6;
    dmft::Rng rng(3);
    Matrix<double> a = random_matrix(rng, n);
    std::vector<double> r = random_vector(rng, n);

    Matrix<double> a_new = a;
    std::vector<double> x_new = r;
    const bool ok_new = dmft::solve_in_place(a_new, x_new);

    std::vector<double> a_legacy_flat = flatten(a);
    std::vector<double> x_legacy = r;
    const bool ok_legacy = legacy_bridge::solve_in_place_double(a_legacy_flat, x_legacy, n);

    expect(ok_new == ok_legacy, "solve_in_place<double> agrees with legacy div_left() on success");
    expect(max_abs_diff(x_new, x_legacy) < 1e-9, "solve_in_place<double> matches legacy div_left()");
}

void test_conjugate_gradient() {
    const int n = 5;
    dmft::Rng rng(4);
    Matrix<double> b = random_matrix(rng, n);
    // A = B^T B + n*I: symmetric positive definite by construction.
    Matrix<double> a(n, n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double s = 0.0;
            for (int k = 0; k < n; ++k) s += b[k][i] * b[k][j];
            a[i][j] = s + (i == j ? static_cast<double>(n) : 0.0);
        }
    std::vector<double> rhs = random_vector(rng, n);

    std::vector<double> x_new;
    const bool ok_new = dmft::conjugate_gradient(a, rhs, x_new);

    std::vector<double> a_flat = flatten(a);
    std::vector<double> x_legacy;
    const bool ok_legacy = legacy_bridge::conjugate_gradient(a_flat, rhs, x_legacy, n);

    expect(ok_new && ok_legacy, "conjugate_gradient converges (new and legacy)");
    expect(max_abs_diff(x_new, x_legacy) < 1e-9, "conjugate_gradient matches legacy div_conj_grad()");
}

void test_jacobi_eigensolver() {
    const int n = 6;
    dmft::Rng rng(5);
    Matrix<double> h = random_symmetric_matrix(rng, n);

    Matrix<double> vecs_new;
    std::vector<double> vals_new;
    dmft::jacobi_eigensolver(h, vecs_new, vals_new);

    std::vector<double> h_flat = flatten(h);
    std::vector<double> vecs_legacy_flat, vals_legacy;
    legacy_bridge::jacobi_eigensolver(h_flat, vecs_legacy_flat, vals_legacy, n);

    expect(max_abs_diff(vals_new, vals_legacy) < 1e-9,
           "jacobi_eigensolver eigenvalues match legacy EigenJacobi()");
    expect(max_abs_diff(flatten(vecs_new), vecs_legacy_flat) < 1e-9,
           "jacobi_eigensolver eigenvectors match legacy EigenJacobi()");
}

void test_block_eigensolver() {
    const int n = 8;
    const int n_blocks = 3;
    std::vector<int> block = {0, 1, 0, 2, 1, 2, 0, 2};

    dmft::Rng rng(6);
    Matrix<double> h(n, n, 0.0);
    // Only fill entries within a block, so the block structure is exact
    // (an off-block entry would make EigenBlock's result plain wrong, not
    // just a different rewrite of the same algorithm).
    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j) {
            if (block[i] != block[j]) continue;
            const double v = 2.0 * rng.uniform01() - 1.0;
            h[i][j] = v;
            h[j][i] = v;
        }

    Matrix<double> vecs_new;
    std::vector<double> vals_new;
    dmft::block_eigensolver(h, vecs_new, vals_new, block, n_blocks);

    std::vector<double> h_flat = flatten(h);
    std::vector<double> vecs_legacy_flat, vals_legacy;
    legacy_bridge::block_eigensolver(h_flat, vecs_legacy_flat, vals_legacy, block, n, n_blocks);

    expect(max_abs_diff(vals_new, vals_legacy) < 1e-9,
           "block_eigensolver eigenvalues match legacy EigenBlock()");
    expect(max_abs_diff(flatten(vecs_new), vecs_legacy_flat) < 1e-9,
           "block_eigensolver eigenvectors match legacy EigenBlock()");
}

void test_leading_eigenvectors_double() {
    const int n = 6;
    // The legacy Eigen_leading has no separate "how many" parameter: its
    // dimension argument doubles as the count, so a fair comparison needs
    // count == n here (leading_eigenvectors's ability to ask for fewer
    // than n is a genuine extension with no legacy equivalent to compare
    // against).
    const int count = n;
    const unsigned int seed = 42;
    dmft::Rng gen(7);
    Matrix<double> h = random_symmetric_matrix(gen, n);

    Matrix<double> vecs_new(count, n);
    std::vector<double> vals_new(count);
    dmft::Rng rng(seed);
    const int found_new = dmft::leading_eigenvectors(h, vecs_new, vals_new, count, rng);

    std::vector<double> h_flat = flatten(h);
    std::vector<double> vecs_legacy_flat, vals_legacy;
    const int found_legacy = legacy_bridge::leading_eigenvectors_double(
        h_flat, vecs_legacy_flat, vals_legacy, n, count, seed, 1e-10);

    expect(found_new == found_legacy, "leading_eigenvectors<double> agrees on count found");
    std::vector<double> vals_new_padded(vals_legacy.size(), 0.0);
    for (int i = 0; i < count; ++i) vals_new_padded[i] = vals_new[i];
    expect(max_abs_diff(vals_new_padded, vals_legacy) < 1e-9,
           "leading_eigenvectors<double> eigenvalues match legacy Eigen_leading()");

    std::vector<double> vecs_new_flat(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < count; ++i)
        for (int j = 0; j < n; ++j) vecs_new_flat[i * n + j] = vecs_new[i][j];
    expect(max_abs_diff(vecs_new_flat, vecs_legacy_flat) < 1e-9,
           "leading_eigenvectors<double> eigenvectors match legacy Eigen_leading()");
}

void test_leading_eigenvectors_complex() {
    const int n = 6;
    // The legacy Eigen_leading has no separate "how many" parameter: its
    // dimension argument doubles as the count, so a fair comparison needs
    // count == n here (leading_eigenvectors's ability to ask for fewer
    // than n is a genuine extension with no legacy equivalent to compare
    // against).
    const int count = n;
    const unsigned int seed = 43;
    dmft::Rng gen(8);
    // Hermitian, so the power-iteration deflation (which conjugates one
    // side, per the original) is operating on a well-posed problem.
    Matrix<Complex> h(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j) {
            if (i == j) {
                h[i][j] = Complex(2.0 * gen.uniform01() - 1.0, 0.0);
            } else {
                const Complex v(2.0 * gen.uniform01() - 1.0, 2.0 * gen.uniform01() - 1.0);
                h[i][j] = v;
                h[j][i] = std::conj(v);
            }
        }

    Matrix<Complex> vecs_new(count, n);
    std::vector<double> vals_new(count);
    dmft::Rng rng(seed);
    const int found_new = dmft::leading_eigenvectors(h, vecs_new, vals_new, count, rng);

    std::vector<Complex> h_flat = flatten(h);
    std::vector<Complex> vecs_legacy_flat;
    std::vector<double> vals_legacy;
    const int found_legacy = legacy_bridge::leading_eigenvectors_complex(
        h_flat, vecs_legacy_flat, vals_legacy, n, count, seed, 1e-10);

    expect(found_new == found_legacy, "leading_eigenvectors<Complex> agrees on count found");
    std::vector<double> vals_new_padded(vals_legacy.size(), 0.0);
    for (int i = 0; i < count; ++i) vals_new_padded[i] = vals_new[i];
    expect(max_abs_diff(vals_new_padded, vals_legacy) < 1e-9,
           "leading_eigenvectors<Complex> eigenvalues match legacy Eigen_leading()");

    std::vector<Complex> vecs_new_flat(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < count; ++i)
        for (int j = 0; j < n; ++j) vecs_new_flat[i * n + j] = vecs_new[i][j];
    expect(max_abs_diff(vecs_new_flat, vecs_legacy_flat) < 1e-9,
           "leading_eigenvectors<Complex> eigenvectors match legacy Eigen_leading()");
}

void test_rng_matches_legacy_int_rnd() {
    // The random-number source itself: Rng(seed) must reproduce int_rnd()
    // started from the same INT_RANDOM, since every test above relies on
    // that to make the two sides comparable at all.
    const unsigned int seed = 999;
    const std::vector<long> legacy_draws = legacy_bridge::raw_draws(seed, 8);

    dmft::Rng rng(seed);
    bool all_match = true;
    for (long expected : legacy_draws)
        if (rng.next_raw() != expected) all_match = false;
    expect(all_match, "Rng(seed).next_raw() matches legacy int_rnd() draw-for-draw");
}

}  // namespace

int main() {
    test_invert();
    test_solve_complex();
    test_solve_double();
    test_conjugate_gradient();
    test_jacobi_eigensolver();
    test_block_eigensolver();
    test_leading_eigenvectors_double();
    test_leading_eigenvectors_complex();
    test_rng_matches_legacy_int_rnd();

    std::printf("\n%s\n", g_failures == 0 ? "ALL TESTS PASSED" : "TESTS FAILED");
    return g_failures == 0 ? 0 : 1;
}
