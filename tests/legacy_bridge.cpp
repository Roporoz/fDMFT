#include "legacy_bridge.hpp"

#include "../legacy/headers.new.h"

namespace legacy_bridge {

std::vector<long> raw_draws(unsigned int seed, int count) {
    INT_RANDOM = static_cast<int>(seed);
    std::vector<long> draws(count);
    for (int i = 0; i < count; ++i) draws[i] = int_rnd();
    return draws;
}

void invert(std::vector<Complex>& a_flat, int n) {
    complex** a = new_complex2(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) a[i][j] = a_flat[i * n + j];

    Inverse(a, n);

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) a_flat[i * n + j] = a[i][j];
    delete_complex2(a, n, n);
}

void solve_in_place_complex(std::vector<Complex>& a_flat, std::vector<Complex>& r, int n) {
    complex** a = new_complex2(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) a[i][j] = a_flat[i * n + j];

    complex* r_raw = new complex[n];
    for (int i = 0; i < n; ++i) r_raw[i] = r[i];

    div_left(a, r_raw, n);

    for (int i = 0; i < n; ++i) r[i] = r_raw[i];
    delete[] r_raw;
    delete_complex2(a, n, n);
}

bool solve_in_place_double(std::vector<double>& a_flat, std::vector<double>& r, int n) {
    double** a = new_double2(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) a[i][j] = a_flat[i * n + j];

    double* r_raw = new double[n];
    for (int i = 0; i < n; ++i) r_raw[i] = r[i];

    bool ok = div_left(a, r_raw, n);

    for (int i = 0; i < n; ++i) r[i] = r_raw[i];
    delete[] r_raw;
    delete_double2(a, n, n);
    return ok;
}

bool conjugate_gradient(const std::vector<double>& a_flat, const std::vector<double>& b,
                         std::vector<double>& x, int n) {
    double** a = new_double2(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) a[i][j] = a_flat[i * n + j];

    double* b_raw = new double[n];
    double* x_raw = new double[n];
    for (int i = 0; i < n; ++i) b_raw[i] = b[i];

    bool ok = div_conj_grad(a, b_raw, x_raw, n);

    x.resize(n);
    for (int i = 0; i < n; ++i) x[i] = x_raw[i];
    delete[] b_raw;
    delete[] x_raw;
    delete_double2(a, n, n);
    return ok;
}

void jacobi_eigensolver(const std::vector<double>& h_flat, std::vector<double>& eigenvectors_flat,
                         std::vector<double>& eigenvalues, int n) {
    double** h = new_double2(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) h[i][j] = h_flat[i * n + j];

    double** a = new_double2(n, n);
    double* e = new double[n];

    EigenJacobi(h, a, e, n);

    eigenvectors_flat.resize(static_cast<std::size_t>(n) * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) eigenvectors_flat[i * n + j] = a[i][j];
    eigenvalues.assign(e, e + n);

    delete[] e;
    delete_double2(a, n, n);
    delete_double2(h, n, n);
}

void block_eigensolver(const std::vector<double>& h_flat, std::vector<double>& eigenvectors_flat,
                        std::vector<double>& eigenvalues, const std::vector<int>& block, int n,
                        int n_blocks) {
    double** h = new_double2(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) h[i][j] = h_flat[i * n + j];

    double** a = new_double2(n, n);
    double* e = new double[n];
    int* block_raw = new int[n];
    for (int i = 0; i < n; ++i) block_raw[i] = block[i];

    EigenBlock(h, a, e, n, block_raw, n_blocks);

    eigenvectors_flat.resize(static_cast<std::size_t>(n) * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) eigenvectors_flat[i * n + j] = a[i][j];
    eigenvalues.assign(e, e + n);

    delete[] block_raw;
    delete[] e;
    delete_double2(a, n, n);
    delete_double2(h, n, n);
}

int leading_eigenvectors_double(const std::vector<double>& h_flat,
                                 std::vector<double>& eigenvectors_flat,
                                 std::vector<double>& eigenvalues, int n, int count,
                                 unsigned int seed, double e_min) {
    double** h = new_double2(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) h[i][j] = h_flat[i * n + j];

    double** a = new_double2(n, n);
    double* e = new double[n];

    INT_RANDOM = static_cast<int>(seed);
    int found = Eigen_leading(h, a, e, count, e_min);

    eigenvectors_flat.assign(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < count; ++i)
        for (int j = 0; j < n; ++j) eigenvectors_flat[i * n + j] = a[i][j];
    eigenvalues.assign(e, e + n);

    delete[] e;
    delete_double2(a, n, n);
    delete_double2(h, n, n);
    return found;
}

int leading_eigenvectors_complex(const std::vector<Complex>& h_flat,
                                  std::vector<Complex>& eigenvectors_flat,
                                  std::vector<double>& eigenvalues, int n, int count,
                                  unsigned int seed, double e_min) {
    complex** h = new_complex2(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) h[i][j] = h_flat[i * n + j];

    complex** a = new_complex2(n, n);
    double* e = new double[n];

    INT_RANDOM = static_cast<int>(seed);
    int found = Eigen_leading(h, a, e, count, e_min);

    eigenvectors_flat.assign(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < count; ++i)
        for (int j = 0; j < n; ++j) eigenvectors_flat[i * n + j] = a[i][j];
    eigenvalues.assign(e, e + n);

    delete[] e;
    delete_complex2(a, n, n);
    delete_complex2(h, n, n);
    return found;
}

}  // namespace legacy_bridge
