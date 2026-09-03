#pragma once

// Plain, macro-free declarations for the legacy numerical routines,
// implemented in legacy_bridge.cpp (which is the only place that includes
// legacy/headers.new.h, with its `#define complex complex<double>` and
// `using namespace std;`). Included from test_linalg.cpp, which links
// against dmft/linalg.hpp's rewrite and must not see any of that.
//
// All matrices are passed as flat, row-major vectors of length n*n so this
// header needs nothing beyond <complex> and <vector>.

#include <complex>
#include <vector>

namespace legacy_bridge {

using Complex = std::complex<double>;

// int_rnd(): reseeds INT_RANDOM=seed, then returns `count` successive raw
// draws. Used to check dmft::Rng reproduces the legacy generator exactly.
std::vector<long> raw_draws(unsigned int seed, int count);

// Inverse(complex**&, int): inverts a_flat (n x n) in place.
void invert(std::vector<Complex>& a_flat, int n);

// div_left(complex**, complex*, int): solves a_flat * x = r for x,
// overwriting r. a_flat is destroyed, matching the original's contract.
void solve_in_place_complex(std::vector<Complex>& a_flat, std::vector<Complex>& r, int n);

// div_left(double**, double*, int): as above, real-valued; returns false on
// a singular pivot.
bool solve_in_place_double(std::vector<double>& a_flat, std::vector<double>& r, int n);

// div_conj_grad(double**, double*, double*, int).
bool conjugate_gradient(const std::vector<double>& a_flat, const std::vector<double>& b,
                         std::vector<double>& x, int n);

// EigenJacobi(double**, double**, double*, int): h_flat (n x n, symmetric)
// -> eigenvectors_flat (n x n, row i is the i-th eigenvector) and
// eigenvalues.
void jacobi_eigensolver(const std::vector<double>& h_flat, std::vector<double>& eigenvectors_flat,
                         std::vector<double>& eigenvalues, int n);

// EigenBlock(double**, double**, double*, int, int*, int).
void block_eigensolver(const std::vector<double>& h_flat, std::vector<double>& eigenvectors_flat,
                        std::vector<double>& eigenvalues, const std::vector<int>& block, int n,
                        int n_blocks);

// Eigen_leading(double**, ...): reseeds the legacy global RNG (INT_RANDOM)
// from `seed` immediately before running, for reproducibility.
int leading_eigenvectors_double(const std::vector<double>& h_flat,
                                 std::vector<double>& eigenvectors_flat,
                                 std::vector<double>& eigenvalues, int n, int count,
                                 unsigned int seed, double e_min);

// Eigen_leading(complex**, ...): same seeding contract as the double
// overload above.
int leading_eigenvectors_complex(const std::vector<Complex>& h_flat,
                                  std::vector<Complex>& eigenvectors_flat,
                                  std::vector<double>& eigenvalues, int n, int count,
                                  unsigned int seed, double e_min);

}  // namespace legacy_bridge
