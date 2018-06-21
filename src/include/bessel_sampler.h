/*
 *
 *    Copyright (c) 2016-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_BESSEL_SAMPLER_H_
#define SRC_INCLUDE_BESSEL_SAMPLER_H_

#include <utility>
#include <vector>

#include "logging.h"
#include "random.h"

namespace smash {

/**
 * The intention of this class is to efficiently sample \f$ (N_1, N_2) \f$
 * from the Skellam distribution \f$ p(N_1,N_2) \sim \mathrm{Poi}(\nu_1)
 * \mathrm{Poi}(\nu_2) \delta(N_1 - N_2 = N)\f$, where \f$\mathrm{Poi}(\nu)\f$
 * denotes Poisson distribution with mean \f$\nu\f$. In other words, this
 * class samples two Poisson numbers with a given mean and a fixed difference.
 * The intended use is to sample the number of baryons and antibaryons, given
 * their means and net baryon number.
 *
 * The distribution of \f$ min(N_1,N_2) \f$ is a so-called Bessel distribution.
 * Denoting \f$ a = \sqrt{\nu_1 \nu_2}\f$, \f$ p(N_{smaller} = k) =
 * \frac{(a/2)^{2k+N}}{I_N(a) k! (N+k)!} \f$. We sample this distribution using
 * the method suggested by Yuan and Kalbfleisch \iref{Yuan2000}: if
 * \f$ m = \frac{1}{2} (\sqrt{a^2 + N^2} - N) > 6\f$, then the distribution is
 * approximated well by a Gaussian, else probabilities are computed explicitely
 * and a table sampling is applied.
 */
class BesselSampler {
 public:
  /**
   * Construct a \ref BesselSampler.
   *
   * \param[in] poisson_mean1 Mean of the first number's Poisson distribution.
   * \param[in] poisson_mean2 Mean of the second number's Poisson distribution.
   * \param[in] fixed_difference Difference between the sampled numbers.
   * \return Constructed sampler.
   */
  BesselSampler(const double poisson_mean1, const double poisson_mean2,
                const int fixed_difference)
      : a_(2.0 * std::sqrt(poisson_mean1 * poisson_mean2)),
        N_(std::abs(fixed_difference)),
        N_is_positive_(fixed_difference >= 0) {
    const auto &log = logger<LogArea::GrandcanThermalizer>();
    assert(poisson_mean1 > 0.0);
    assert(poisson_mean2 > 0.0);
    log.debug("Bessel sampler", ": Poisson mean N1 = ", poisson_mean1,
              ", Poisson mean N2 = ", poisson_mean2, ", N1 - N2 fixed to ",
              fixed_difference);
    m_ = 0.5 * (std::sqrt(a_ * a_ + N_ * N_) - N_);
    if (m_ >= m_switch_method_) {
      mu_ = 0.5 * a_ * r(N_, a_);
      const double mean_sqr = mu_ * (1.0 + 0.5 * a_ * r(N_ + 1, a_));
      sigma_ = std::sqrt(mean_sqr - mu_ * mu_);
      log.debug("m = ", m_, " -> using gaussian sampling with mean = ", mu_,
                ", sigma = ", sigma_);
    } else {
      log.debug("m = ", m_, " -> using direct sampling method");
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
        log.debug("Probability (", i, ") = ", p);
        i++;
      }
      dist_.reset_weights(probabilities);
    }
  }

  /**
   * Sample two numbers from given Poissonians with a fixed difference.
   *
   * \return Pair of first and second sampled number.
   */
  std::pair<int, int> sample() {
    const int N_smaller = (m_ >= m_switch_method_)
                              ? std::round(random::normal(mu_, sigma_))
                              : dist_();
    return N_is_positive_ ? std::make_pair(N_smaller + N_, N_smaller)
                          : std::make_pair(N_smaller, N_smaller + N_);
  }

 private:
  /**
   * Compute the ratio of two Bessel functions
   * r(n,a) = bessel_I(n+1,a)/bessel_I(n,a) using the continued fraction
   * representation (see \iref{Yuan2000}).
   *
   * \param[in] n First Bessel parameter.
   * \param[in] a Second Bessel parameter.
   * \return Ratio bessel_I(n+1,a)/bessel_I(n,a).
   */
  static double r(int n, double a) {
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
  /// Vector to store tabulated values of probabilities for small m case (m <6).
  random::discrete_dist<double> dist_;

  /// Mode of the Bessel function, see \iref{Yuan2000} for details.
  double m_;

  /// Second parameter of Bessel distribution, see \iref{Yuan2000} for details.
  const double a_;

  /// First parameter of Bessel distribution (= \f$ \nu \f$ in \iref{Yuan2000}).
  const int N_;

  /// Boolean variable to verify that N > 0.
  const bool N_is_positive_;

  /**
   * Switching mode to normal approximation.
   * \note Normal approximation of Bessel functions is possible for modes >= 6.
   * See \iref{Yuan2000} for details.
   */
  static constexpr double m_switch_method_ = 6.0;

  /// Probabilities smaller than negligibly_probability are neglected.
  static constexpr double negligible_probability_ = 1.e-12;

  /// Mean of the Bessel distribution.
  double mu_;

  /// Standard deviation of the Bessel distribution.
  double sigma_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_BESSEL_SAMPLER_H_
