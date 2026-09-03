#pragma once

#include <cmath>
#include <complex>

// Replaces legacy headers.new.h's `#define complex complex<double>` (which
// shadowed std::complex itself and broke any header included afterwards
// that used "complex" as an identifier) and its free-standing `I`/`Pi`
// globals.
namespace dmft {

using Complex = std::complex<double>;

inline const double pi = 4.0 * std::atan(1.0);
inline const Complex I{0.0, 1.0};

template <typename T>
int sgn(T x) {
    if (x > T{0}) return 1;
    if (x < T{0}) return -1;
    return 0;
}

template <typename T>
T sqr(T x) {
    return x * x;
}

inline Complex sqr(const Complex& x) { return x * x; }

// |x|^2, matching the original norm2() (avoids the sqrt in std::abs).
inline double norm2(const Complex& x) {
    return x.real() * x.real() + x.imag() * x.imag();
}

}  // namespace dmft
