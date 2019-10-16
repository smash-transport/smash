/*
 *
 *    Copyright (c) 2014-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/random.h"
#include <random>
#include "smash/logging.h"

namespace smash {
static constexpr int LGrandcanThermalizer = LogArea::GrandcanThermalizer::id;
/*thread_local (see #3075)*/ random::Engine random::engine;

int64_t random::generate_63bit_seed() {
  std::random_device rd;
  static_assert(std::is_same<decltype(rd()), uint32_t>::value,
                "random_device is assumed to generate uint32_t");
  uint64_t unsigned_seed =
      (static_cast<uint64_t>(rd()) << 32) | static_cast<uint64_t>(rd());
  // Discard the highest bit to make sure it fits into a positive int64_t
  const int64_t seed = static_cast<int64_t>(unsigned_seed >> 1);
  return seed;
}

random::BesselSampler::BesselSampler(const double poisson_mean1,
                                     const double poisson_mean2,
                                     const int fixed_difference)
    : a_(2.0 * std::sqrt(poisson_mean1 * poisson_mean2)),
      N_(std::abs(fixed_difference)),
      N_is_positive_(fixed_difference >= 0) {
  assert(poisson_mean1 >= 0.0);
  assert(poisson_mean2 >= 0.0);
  logg[LGrandcanThermalizer].debug("Bessel sampler",
                                   ": Poisson mean N1 = ", poisson_mean1,
                                   ", Poisson mean N2 = ", poisson_mean2,
                                   ", N1 - N2 fixed to ", fixed_difference);
  m_ = 0.5 * (std::sqrt(a_ * a_ + N_ * N_) - N_);
  if (m_ >= m_switch_method_) {
    mu_ = 0.5 * a_ * r_(N_, a_);
    const double mean_sqr = mu_ * (1.0 + 0.5 * a_ * r_(N_ + 1, a_));
    sigma_ = std::sqrt(mean_sqr - mu_ * mu_);
    logg[LGrandcanThermalizer].debug(
        "m = ", m_, " -> using gaussian sampling with mean = ", mu_,
        ", sigma = ", sigma_);
  } else {
    logg[LGrandcanThermalizer].debug("m = ", m_,
                                     " -> using direct sampling method");
    std::vector<double> probabilities;
    double wi = 1.0, sum = 0.0;
    int i = 0;
    do {
      sum += wi;
      probabilities.push_back(wi);
      wi *= 0.25 * a_ * a_ / (i + 1) / (N_ + i + 1);
      i++;
    } while (wi > negligible_probability_);
    i = 0;
    for (double p : probabilities) {
      p /= sum;
      logg[LGrandcanThermalizer].debug("Probability (", i, ") = ", p);
      i++;
    }
    dist_.reset_weights(probabilities);
  }
}

std::pair<int, int> random::BesselSampler::sample() {
  const int N_smaller = (m_ >= m_switch_method_)
                            ? std::round(random::normal(mu_, sigma_))
                            : dist_();
  return N_is_positive_ ? std::make_pair(N_smaller + N_, N_smaller)
                        : std::make_pair(N_smaller, N_smaller + N_);
}

double random::BesselSampler::r_(int n, double a) {
  const double a_inv = 1.0 / a;
  double res = 0.0;
  // |x - continued fraction of order n| < 2^(-n+1), see the book
  // "Continued fractions" by Khinchin. For 10^-16 = ~2^-50 precision
  // 50 iterations should be sufficient. However, I found that for some
  // numerical reason at least 100 terms are needed.
  int i = 200;
  for (; i > 0; i--) {
    res = 1.0 / (a_inv * 2 * (n + i) + res);
  }
  // Check the known property of r(n,a) function, see iref{Yuan2000}.
  assert(a / (std::sqrt(a * a + (n + 1) * (n + 1)) + n + 1) <= res);
  assert(res <= a / (std::sqrt(a * a + n * n) + n));
  return res;
}

}  // namespace smash
