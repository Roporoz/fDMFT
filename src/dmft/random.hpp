#pragma once

#include <complex>
#include <ctime>
#include <stdexcept>

#include "core.hpp"

// Replaces legacy headers.new.h's int_rnd()/rnd()/rnd_gauss2(), which lived
// as free functions mutating a single process-wide global (INT_RANDOM). That
// made every draw depend on global call order and made the whole program
// non-reproducible between runs (the seed came from time(NULL)) and
// impossible to run two independent streams side by side.
//
// Rng is the same Lehmer/Park-Miller generator (multiplier 16807, modulus
// 2^31-1, Schrage's method to avoid overflow) wrapped as an explicit,
// seedable, reentrant object: two Rng instances never interfere, and a
// fixed seed reproduces a run exactly.
namespace dmft {

class Rng {
public:
    // Matches the original module-load-time behavior: seed from the clock,
    // then discard 4 draws ("randomising..." in the original comment) so
    // consecutive processes started in the same second still diverge.
    Rng() : Rng(static_cast<unsigned int>(std::time(nullptr))) { warm_up(); }

    // Deterministic seeding for reproducible runs and tests. No warm-up is
    // applied here (unlike the time-seeded constructor above) so that a
    // given seed always maps to the same first draw.
    explicit Rng(unsigned int seed) : state_(seed == 0 ? 1 : static_cast<long>(seed % modulus)) {}

    // Uniform integer in [1, 2^31 - 2]. Equivalent to the original int_rnd().
    long next_raw() {
        long k = state_ / 127773;
        state_ = 16807 * (state_ - k * 127773) - 2836 * k;
        if (state_ < 0) state_ += modulus;
        return state_;
    }

    // Uniform double in (0, 1). Equivalent to the original rnd().
    double uniform01() { return static_cast<double>(next_raw()) / static_cast<double>(modulus); }

    // Uniform integer in [0, k). Equivalent to the original rnd(int).
    int uniform_int(int k) {
        if (k <= 0) throw std::invalid_argument("Rng::uniform_int: k must be positive");
        const long d = modulus % k;
        const long d1 = modulus - d;
        long r;
        do {
            r = next_raw();
        } while (r >= d1);
        return static_cast<int>(r / (d1 / k));
    }

    // Complex draw whose squared magnitude is Exp(1)-distributed with a
    // uniform phase (a scaled 2D Gaussian-like proposal, not literally
    // Box-Muller). Equivalent to the original rnd_gauss2(). Used only as a
    // jump proposal in stochastic bath fitting / initial-vector seeding, so
    // the exact statistical shape doesn't affect physical results.
    Complex gaussian2d() {
        const double R = std::sqrt(-std::log(uniform01()));
        const double phi = 2.0 * pi * uniform01();
        return R * std::exp(I * phi);
    }

private:
    static constexpr long modulus = 2147483647;  // 2^31 - 1

    void warm_up() {
        for (int i = 0; i < 4; ++i) next_raw();
    }

    long state_;
};

}  // namespace dmft
