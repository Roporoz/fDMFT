#include "linalg.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace dmft {

void invert(Matrix<Complex>& a) {
    const int size = a.n_rows();
    Matrix<Complex> r = Matrix<Complex>::identity(size);

    for (int i = 0; i < size; ++i) {
        double fmax = -1.0;
        int jmax = i;
        for (int j = i; j < size; ++j) {
            if (fmax < norm2(a[j][i])) {
                jmax = j;
                fmax = norm2(a[j][i]);
            }
        }
        std::swap(a[i], a[jmax]);
        std::swap(r[i], r[jmax]);
        if (norm2(a[i][i]) <= 0.0) throw std::runtime_error("invert: matrix is singular");

#pragma omp parallel for
        for (int j = 0; j < i; ++j) {
            const Complex f = a[j][i] / a[i][i];
            a[j][i] = 0.0;
            for (int l = i + 1; l < size; ++l) a[j][l] -= a[i][l] * f;
            for (int l = 0; l < size; ++l) r[j][l] -= r[i][l] * f;
        }
        {
            const Complex f = 1.0 / a[i][i];
            for (int l = i; l < size; ++l) a[i][l] *= f;
            for (int l = 0; l < size; ++l) r[i][l] *= f;
        }
#pragma omp parallel for
        for (int j = i + 1; j < size; ++j) {
            const Complex f = a[j][i] / a[i][i];
            a[j][i] = 0.0;
            for (int l = i + 1; l < size; ++l) a[j][l] -= a[i][l] * f;
            for (int l = 0; l < size; ++l) r[j][l] -= r[i][l] * f;
        }
    }

    a = std::move(r);
}

bool solve_in_place(Matrix<Complex>& a, std::vector<Complex>& r) {
    const int size = a.n_rows();
    for (int i = 0; i < size; ++i) {
        double fmax = -1.0;
        int jmax = i;
        for (int j = i; j < size; ++j) {
            if (fmax < norm2(a[j][i])) {
                jmax = j;
                fmax = norm2(a[j][i]);
            }
        }
        std::swap(a[i], a[jmax]);
        std::swap(r[i], r[jmax]);
        if (norm2(a[i][i]) <= 0.0) return false;

#pragma omp parallel for
        for (int j = 0; j < i; ++j) {
            const Complex f = a[j][i] / a[i][i];
            a[j][i] = 0.0;
            for (int l = i + 1; l < size; ++l) a[j][l] -= a[i][l] * f;
            r[j] -= r[i] * f;
        }
        {
            const Complex f = 1.0 / a[i][i];
            for (int l = i; l < size; ++l) a[i][l] *= f;
            r[i] *= f;
        }
#pragma omp parallel for
        for (int j = i + 1; j < size; ++j) {
            const Complex f = a[j][i] / a[i][i];
            a[j][i] = 0.0;
            for (int l = i + 1; l < size; ++l) a[j][l] -= a[i][l] * f;
            r[j] -= r[i] * f;
        }
    }
    return true;
}

bool solve_in_place(Matrix<double>& a, std::vector<double>& r) {
    const int size = a.n_rows();
    for (int i = 0; i < size; ++i) {
        double fmax = -1.0;
        int jmax = i;
        for (int j = i; j < size; ++j) {
            if (fmax < std::abs(a[j][i])) {
                jmax = j;
                fmax = std::abs(a[j][i]);
            }
        }
        std::swap(a[i], a[jmax]);
        std::swap(r[i], r[jmax]);
        if (std::abs(a[i][i]) <= 0.0) return false;

#pragma omp parallel for
        for (int j = 0; j < i; ++j) {
            const double f = a[j][i] / a[i][i];
            a[j][i] = 0.0;
            for (int l = i + 1; l < size; ++l) a[j][l] -= a[i][l] * f;
            r[j] -= r[i] * f;
        }
        {
            const double f = 1.0 / a[i][i];
            for (int l = i; l < size; ++l) a[i][l] *= f;
            r[i] *= f;
        }
#pragma omp parallel for
        for (int j = i + 1; j < size; ++j) {
            const double f = a[j][i] / a[i][i];
            a[j][i] = 0.0;
            for (int l = i + 1; l < size; ++l) a[j][l] -= a[i][l] * f;
            r[j] -= r[i] * f;
        }
    }
    return true;
}

bool conjugate_gradient(const Matrix<double>& A, const std::vector<double>& b,
                         std::vector<double>& x, int max_iterations, double tol) {
    const int size = A.n_rows();
    if (max_iterations < 0) max_iterations = 5 * size;

    std::vector<double> p(size), r(size), Ap(size);
    x.assign(size, 0.0);
    for (int j = 0; j < size; ++j) {
        r[j] = b[j];
        p[j] = b[j];
    }

    for (int i = 0; i < max_iterations; ++i) {
        double r2 = 0.0;
        for (int j = 0; j < size; ++j) r2 += r[j] * r[j];

        for (int j = 0; j < size; ++j) {
            Ap[j] = 0.0;
            for (int k = 0; k < size; ++k) Ap[j] += A[j][k] * p[k];
        }

        double p2 = 0.0;
        for (int j = 0; j < size; ++j) p2 += p[j] * Ap[j];
        const double alpha = r2 / p2;
        for (int j = 0; j < size; ++j) x[j] += alpha * p[j];

        for (int j = 0; j < size; ++j) r[j] -= alpha * Ap[j];
        double r2_new = 0.0;
        for (int j = 0; j < size; ++j) r2_new += r[j] * r[j];
        if (r2_new < tol) return true;

        const double beta = r2_new / r2;
        for (int j = 0; j < size; ++j) p[j] = r[j] + beta * p[j];
    }
    return false;
}

void jacobi_eigensolver(const Matrix<double>& h, Matrix<double>& eigenvectors,
                         std::vector<double>& eigenvalues, double tol, int max_sweeps) {
    const int n = h.n_rows();
    Matrix<double> a(n, n);
    Matrix<double> v = Matrix<double>::identity(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) a[i][j] = h[i][j];

    bool changed;
    do {
        --max_sweeps;
        changed = false;
        for (int pp = n - 1; pp >= 1; --pp) {
            for (int q = 0; q < n - pp; ++q) {
                const int p = pp + q;
                if (std::abs(a[p][q]) <= tol) continue;
                changed = true;

                const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                double t = 1.0 / (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                if (theta < 0.0) t = -t;
                const double c = 1.0 / std::sqrt(t * t + 1.0);
                const double s = t * c;
                const double tau = s / (1.0 + c);

                for (int r = 0; r < n; ++r) {
                    if (r != p && r != q) {
                        const double a1 = a[r][p] - s * (a[r][q] + tau * a[r][p]);
                        const double a2 = a[r][q] + s * (a[r][p] - tau * a[r][q]);
                        a[r][p] = a1;
                        a[p][r] = a1;
                        a[r][q] = a2;
                        a[q][r] = a2;
                    }
                    const double v1 = c * v[p][r] - s * v[q][r];
                    const double v2 = s * v[p][r] + c * v[q][r];
                    v[p][r] = v1;
                    v[q][r] = v2;
                }
                const double ap = a[p][p] - t * a[p][q];
                const double aq = a[q][q] + t * a[p][q];
                a[p][p] = ap;
                a[q][q] = aq;
                a[p][q] = 0.0;
                a[q][p] = 0.0;
            }
        }
    } while (changed && max_sweeps > 0);

    eigenvalues.assign(n, 0.0);
    for (int i = 0; i < n; ++i) eigenvalues[i] = a[i][i];

    eigenvectors = std::move(v);
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int j = 0; j < n; ++j) s += eigenvectors[i][j] * eigenvectors[i][j];
        s = std::sqrt(s);
        if (std::abs(s) > tol * tol) {
            for (int j = 0; j < n; ++j) eigenvectors[i][j] /= s;
        } else {
            const double uniform = 1.0 / std::sqrt(static_cast<double>(n));
            for (int j = 0; j < n; ++j) eigenvectors[i][j] = uniform;
        }
    }
}

void block_eigensolver(const Matrix<double>& h, Matrix<double>& eigenvectors,
                        std::vector<double>& eigenvalues, const std::vector<int>& block,
                        int n_blocks, double tol, int max_sweeps) {
    const int n = h.n_rows();
    eigenvectors = Matrix<double>(n, n, 0.0);
    eigenvalues.assign(n, 0.0);

    std::vector<int> members(n);
    for (int b = 0; b < n_blocks; ++b) {
        int nb = 0;
        for (int j = 0; j < n; ++j)
            if (block[j] == b) members[nb++] = j;
        if (nb == 0) continue;

        Matrix<double> hb(nb, nb);
        for (int j1 = 0; j1 < nb; ++j1)
            for (int j2 = 0; j2 < nb; ++j2) hb[j1][j2] = h[members[j1]][members[j2]];

        Matrix<double> ab(nb, nb);
        std::vector<double> eb(nb);
        jacobi_eigensolver(hb, ab, eb, tol, max_sweeps);

        for (int j1 = 0; j1 < nb; ++j1)
            for (int j2 = 0; j2 < nb; ++j2) eigenvectors[members[j1]][members[j2]] = ab[j1][j2];
        for (int j = 0; j < nb; ++j) eigenvalues[members[j]] = eb[j];
    }
}

int leading_eigenvectors(const Matrix<double>& h, Matrix<double>& eigenvectors,
                          std::vector<double>& eigenvalues, int count, Rng& rng, double e_min) {
    const int n = h.n_rows();
    std::vector<double> psi(n), psi_next(n);

    for (int n0 = 0; n0 < count; ++n0) {
        for (int j = 0; j < n; ++j) psi[j] = rng.uniform01() - 0.5;
        double s = 0.0;

        for (int iter = 0; iter < 10000; ++iter) {
            for (int j = 0; j < n; ++j) {
                psi_next[j] = 0.0;
                for (int j1 = 0; j1 < n; ++j1) psi_next[j] += h[j][j1] * psi[j1];
            }
            for (int j = 0; j < n0; ++j) {
                double x = 0.0;
                for (int j1 = 0; j1 < n; ++j1) x += eigenvectors[j][j1] * psi_next[j1];
                for (int j1 = 0; j1 < n; ++j1) psi_next[j1] -= x * eigenvectors[j][j1];
            }
            double norm = 0.0;
            for (int j = 0; j < n; ++j) norm += psi_next[j] * psi_next[j];
            norm = std::sqrt(norm);
            for (int j = 0; j < n; ++j) psi_next[j] /= norm;

            double r = 0.0, r_ = 0.0;
            for (int j = 0; j < n; ++j) {
                r += sqr(psi[j] - psi_next[j]);
                r_ += sqr(psi[j] + psi_next[j]);
                psi[j] = psi_next[j];
            }
            s = norm;
            if (r_ < r) s = -s;
            if (r < 1e-20 || r_ < 1e-20) break;
        }

        for (int j = 0; j < n; ++j) eigenvectors[n0][j] = psi[j];
        eigenvalues[n0] = s;
        if (std::abs(s) < e_min) return n0 + 1;
    }
    return count;
}

int leading_eigenvectors(const Matrix<Complex>& h, Matrix<Complex>& eigenvectors,
                          std::vector<double>& eigenvalues, int count, Rng& rng, double e_min) {
    const int n = h.n_rows();
    std::vector<Complex> psi(n), psi_next(n);

    for (int n0 = 0; n0 < count; ++n0) {
        for (int j = 0; j < n; ++j) psi[j] = Complex(rng.uniform01() - 0.5, rng.uniform01() - 0.5);
        double s = 0.0;

        for (int iter = 0; iter < 10000; ++iter) {
            for (int j = 0; j < n; ++j) {
                psi_next[j] = 0.0;
                for (int j1 = 0; j1 < n; ++j1) psi_next[j] += h[j][j1] * psi[j1];
            }
            for (int j = 0; j < n0; ++j) {
                Complex x = 0.0;
                for (int j1 = 0; j1 < n; ++j1) x += std::conj(eigenvectors[j][j1]) * psi_next[j1];
                for (int j1 = 0; j1 < n; ++j1) psi_next[j1] -= x * eigenvectors[j][j1];
            }
            double norm = 0.0;
            for (int j = 0; j < n; ++j) norm += norm2(psi_next[j]);
            norm = std::sqrt(norm);
            for (int j = 0; j < n; ++j) psi_next[j] /= norm;

            double r = 0.0, r_ = 0.0;
            for (int j = 0; j < n; ++j) {
                r += norm2(psi[j] - psi_next[j]);
                r_ += norm2(psi[j] + psi_next[j]);
                psi[j] = psi_next[j];
            }
            s = norm;
            if (r_ < r) s = -s;
            if (r < 1e-18 || r_ < 1e-18) break;
        }

        for (int j = 0; j < n; ++j) eigenvectors[n0][j] = psi[j];
        eigenvalues[n0] = s;
        if (std::abs(s) < e_min) return n0 + 1;
    }
    return count;
}

}  // namespace dmft
