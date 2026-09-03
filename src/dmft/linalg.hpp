#pragma once

#include <vector>

#include "core.hpp"
#include "random.hpp"

// Replaces legacy headers.new.h's T** allocators (new_double2/new_complex2,
// and new_ldouble2 -- a byte-for-byte duplicate of new_double2 under a
// misleading name, itself a symptom of the copy-paste style this rewrite is
// meant to move away from) together with the free functions that operated
// on them: Inverse, div_left (x2), div_conj_grad, EigenJacobi, EigenBlock,
// Eigen_leading (x2).
//
// Every function below reproduces the original's algorithm and operation
// order exactly (same pivoting rule, same elimination order, same rotation
// formulas), so that for identical inputs it reproduces identical floating
// point results -- see tests/test_linalg.cpp, which links this against the
// original implementations and compares them directly. What changed is the
// interface: RAII storage instead of manual new/delete, explicit Rng
// arguments instead of a hidden global generator, and return values/thrown
// exceptions instead of printing to cout and continuing in an undefined
// state.
namespace dmft {

// A dense, owning matrix. Backed by a vector of row-vectors (mirroring the
// original's T** layout: an array of independently-owned rows) rather than
// one flat buffer, because several algorithms below (partial-pivoting
// elimination, Jacobi's row/column rotations) swap or scale whole rows --
// operations that are O(1) pointer swaps on rows-of-vectors, exactly as
// they were O(1) pointer swaps on the original's T** rows.
template <typename T>
class Matrix {
public:
    Matrix() = default;
    Matrix(int rows, int cols) : rows_(rows, std::vector<T>(static_cast<std::size_t>(cols))) {}
    Matrix(int rows, int cols, const T& fill)
        : rows_(rows, std::vector<T>(static_cast<std::size_t>(cols), fill)) {}

    int n_rows() const { return static_cast<int>(rows_.size()); }
    int n_cols() const { return rows_.empty() ? 0 : static_cast<int>(rows_[0].size()); }

    std::vector<T>& operator[](int i) { return rows_[static_cast<std::size_t>(i)]; }
    const std::vector<T>& operator[](int i) const { return rows_[static_cast<std::size_t>(i)]; }

    static Matrix identity(int n) {
        Matrix m(n, n, T{});
        for (int i = 0; i < n; ++i) m[i][i] = T{1};
        return m;
    }

private:
    std::vector<std::vector<T>> rows_;
};

// Inverts `a` in place (Gauss-Jordan elimination, partial pivoting by
// |a[j][i]|^2 on column i). Throws std::runtime_error if `a` is singular
// (a pivot column is exactly zero) -- the original printed "LInverse!" to
// cout and returned, silently leaving `a` half-eliminated.
// Equivalent to the original Inverse(complex**&, int).
void invert(Matrix<Complex>& a);

// Solves a*x = r for x, overwriting r with the solution. Destroys `a` (left
// row-reduced, not restored) -- same contract as the original. Returns
// false if `a` is singular, leaving r's contents undefined, instead of the
// original complex overload's silent "print and continue" (the original's
// real-valued overload already returned bool for this; this rewrite gives
// both the same contract).
// Equivalent to the original div_left(complex**, complex*, int) and
// div_left(double**, double*, int).
bool solve_in_place(Matrix<Complex>& a, std::vector<Complex>& r);
bool solve_in_place(Matrix<double>& a, std::vector<double>& r);

// Solves the symmetric positive-definite system A*x = b by conjugate
// gradient, for at most max_iterations steps (default 5*size, matching the
// original's fixed budget) or until the residual norm^2 drops below `tol`.
// Returns whether it converged.
// Equivalent to the original div_conj_grad.
bool conjugate_gradient(const Matrix<double>& A, const std::vector<double>& b,
                         std::vector<double>& x, int max_iterations = -1, double tol = 1e-20);

// Diagonalizes the symmetric matrix h via the cyclic Jacobi eigenvalue
// algorithm. On return, eigenvalues[i] is the i-th eigenvalue and
// eigenvectors[i] is its (normalized) eigenvector, in the original's
// (unsorted) order.
// Equivalent to the original EigenJacobi.
void jacobi_eigensolver(const Matrix<double>& h, Matrix<double>& eigenvectors,
                         std::vector<double>& eigenvalues, double tol = 1e-10,
                         int max_sweeps = 100);

// Diagonalizes h assuming it is block-diagonal under the partition given by
// `block` (block[j] is the block index, in [0, n_blocks), of basis state
// j) -- i.e. diagonalizes each block independently via jacobi_eigensolver
// and scatters the results back. Cheaper than a full diagonalization when
// h actually respects the given block structure (e.g. conserved quantum
// numbers), and wrong if it doesn't.
// Equivalent to the original EigenBlock.
void block_eigensolver(const Matrix<double>& h, Matrix<double>& eigenvectors,
                        std::vector<double>& eigenvalues, const std::vector<int>& block,
                        int n_blocks, double tol = 1e-10, int max_sweeps = 100);

// Finds up to `count` leading eigenvectors of h by power iteration with
// deflation against previously found vectors, stopping early (and
// returning fewer than `count`) once a found eigenvalue's magnitude drops
// below e_min. `eigenvectors` and `eigenvalues` must already be sized for
// up to `count` results. Draws random restart vectors from `rng`, where the
// original drew from a hidden global generator.
// Equivalent to the original Eigen_leading(double**, ...) and
// Eigen_leading(complex**, ...).
int leading_eigenvectors(const Matrix<double>& h, Matrix<double>& eigenvectors,
                          std::vector<double>& eigenvalues, int count, Rng& rng,
                          double e_min = 1e-10);
int leading_eigenvectors(const Matrix<Complex>& h, Matrix<Complex>& eigenvectors,
                          std::vector<double>& eigenvalues, int count, Rng& rng,
                          double e_min = 1e-10);

}  // namespace dmft
