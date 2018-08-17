/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_RANDOM_H_
#define SRC_INCLUDE_RANDOM_H_

#include <cassert>
#include <limits>
#include <random>
#include <utility>
#include <vector>

namespace smash {

/** Namespace random provides functions for random Number Generation.
 */

namespace random {

/// The random number engine used is the Mersenne Twister.
using Engine = std::mt19937_64;

/// The engine that is used commonly by all distributions.
extern /*thread_local (see #3075)*/ Engine engine;

/** Provides uniform random numbers on a fixed interval.
 *
 * objects of uniform_dist can be used to provide a large number of
 * random numbers in the same interval. Example:
 *
 * \code
 *   using namespace random;
 *   double sum = 0.0;
 *   auto uniform_0_to_3 = uniform_dist(0., 3.);
 *   for (MANY_TIMES) {
 *     sum += uniform_0_to_3();
 *   }
 * \endcode
 *
 * The random number engine is completely hidden inside the object.
 */
template <typename T>
class uniform_dist {
 public:
  /**
   * Creates the object and fixes the interval.
   *
   * \param min Lower bound of interval.
   * \param max Upper bound of interval.
   * */
  uniform_dist(T min, T max) : distribution(min, max) {}
  /** \returns A random number in the interval. */
  T operator()() { return distribution(engine); }

 private:
  /** The distribution object that is being used. */
  std::uniform_real_distribution<T> distribution;
};

/** Sets the seed of the random number engine. */
template <typename T>
void set_seed(T &&seed) {
  static_assert(std::is_same<Engine::result_type, uint64_t>::value,
                "experiment.cc needs the seed to be 64 bits");
  engine.seed(std::forward<T>(seed));
}

/// Advance the engine's state and return the generated value.
inline Engine::result_type advance() { return engine(); }

/**
 * \returns A uniformly distributed random real number \f$\chi \in [{\rm
 * min}, {\rm max})\f$
 *
 * \param min Minimal sampled value.
 * \param max Maximal sampled value.
 */
template <typename T>
T uniform(T min, T max) {
  return std::uniform_real_distribution<T>(min, max)(engine);
}

/**
 * \return A uniformly distributed random integer number \f$\chi \in [{\rm
 * min}, {\rm max})\f$
 *
 * \param min Minimal sampled value.
 * \param max Maximal sampled value.
 */
template <typename T>
T uniform_int(T min, T max) {
  return std::uniform_int_distribution<T>(min, max)(engine);
}

/**
 * \return a uniformly distributed random number \f$\chi \in [0,1)\f$.
 *
 * Note that the popular implementations in GCC and clang may return 1:
 *
 * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=64351
 * https://llvm.org/bugs/show_bug.cgi?id=18767
 */
template <typename T = double>
T canonical() {
  return std::generate_canonical<T, std::numeric_limits<double>::digits>(
      engine);
}

/**
 * \return A uniformly distributed random number \f$\chi \in (0,1]\f$.
 */
template <typename T = double>
T canonical_nonzero() {
  // use 'nextafter' to generate a value that is guaranteed to be larger than 0
  return std::nextafter(
      std::generate_canonical<T, std::numeric_limits<double>::digits>(engine),
      T(1));
}

/**
 * \return A uniform_dist object.
 * \param min Lower bound of interval.
 * \param max Upper bound of interval.
 * */
template <typename T>
uniform_dist<T> make_uniform_distribution(T min, T max) {
  return uniform_dist<T>(min, max);
}

/**
 * Draws an exponentially distributed random number.
 *
 * Probability for a given return value \f$\chi\f$ is \f$p(\chi) =
 * \Theta(\chi) \cdot \exp(-t)\f$.
 *
 * \param lambda Rate parameter.
 * \return Sampled random number.
 */
template <typename T = double>
T exponential(T lambda) {
  // We are not using std::exponential_distribution because of a bug in the
  // implementations by clang and gcc.
  return -std::log(canonical_nonzero()) / lambda;
}

/**
 * Draws a random number x from an exponential distribution exp(A*x), where A is
 * assumed to be positive, and x is typically negative.
 * The result x is restricted to lie between x1 and x2 (with x2 < x <= x1).
 *
 * \param A Positive shape parameter.
 * \param x1 Maximal sampled value.
 * \param x2 Minimal sampled value.
 * \return Sampled random number.
 */
template <typename T = double>
T expo(T A, T x1, T x2) {
  const T a1 = A * x1, a2 = A * x2;
  const T a_min = std::log(std::numeric_limits<T>::min());
#ifndef NDEBUG
  assert(A > T(0.) && x1 >= x2 && a1 > a_min);
#endif
  const T r1 = std::exp(a1);
  const T r2 = a2 > a_min ? std::exp(a2) : T(0.);  // prevent underflow
  T x;
  do {
    /* sample repeatedly until x is in the requested range
     * (it can get outside due to numerical errors, see issue #2959) */
    x = std::log(uniform(r1, r2)) / A;
  } while (!(x <= x1 && x > x2));
  return x;
}

/**
 * Signum function.
 *
 * \param val The input value.
 * \return The sign of the input value.
 */
template <typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

/**
 * Draws a random number according to a power-law distribution ~ x^n.
 *
 * \param n Exponent in power law (arbitrary real number).
 * \param xMin Minimum value.
 * \param xMax Maximum value.
 * \return random number between xMin and xMax.
 */
template <typename T = double>
T power(T n, T xMin, T xMax) {
  const T n1 = n + 1;
  if (std::abs(n1) < 1E-3) {
    return xMin * std::pow(xMax / xMin, canonical());
  } else if (xMin > 0. && xMax > 0.) {
    return std::pow(uniform(std::pow(xMin, n1), std::pow(xMax, n1)), 1. / n1);
  } else {
    return sgn(xMin) * std::pow(uniform(std::pow(std::abs(xMin), n1),
                                        std::pow(std::abs(xMax), n1)),
                                1. / n1);
  }
}

/**
 * Returns a Poisson distributed random number.
 *
 * Probability for a given return value \f$\chi\f$ is \f$p(\chi) =
 * \chi^i/i! \cdot \exp(-\chi)\f$
 *
 * \param lam Mean value of the distribution.
 * \return Sampled random number.
 */
template <typename T>
int poisson(const T &lam) {
  return std::poisson_distribution<int>(lam)(engine);
}

/**
 * Returns a binomially distributed random number.
 *
 * \param N Number of trials.
 * \param p Probability of a trial generating true.
 * \return Sampled random number.
 */
template <typename T>
int binomial(const int N, const T &p) {
  return std::binomial_distribution<int>(N, p)(engine);
}

/**
 * Returns a random number drawn from a normal distribution.
 *
 * \param mean Mean value of the distribution.
 * \param sigma Standard deviation of the distribution.
 * \return Sampled random number.
 */
template <typename T>
double normal(const T &mean, const T &sigma) {
  return std::normal_distribution<double>(mean, sigma)(engine);
}

/**
 * Discrete distribution with weight given by probability vector.
 */
template <typename T>
class discrete_dist {
 public:
  /** Default discrete distribution.
   *
   * Always draws 0.
   */
  discrete_dist() : distribution({1.0}) {}

  /** Construct from probability vector.
   * \param plist Vector with probabilities such that P(i) = vec[i]
   */
  explicit discrete_dist(const std::vector<T> &plist)
      : distribution(plist.begin(), plist.end()) {}

  /** Construct from probability list.
   * \param l Initializer list with probabilities such that P(i) = l[i]
   */
  explicit discrete_dist(std::initializer_list<T> l) : distribution(l) {}

  /** Reset the discrete distribution from a new probability list.
   * \param plist Vector with probabilities such that P(i) = vec[i]
   */
  void reset_weights(const std::vector<T> &plist) {
    distribution = std::discrete_distribution<>(plist.begin(), plist.end());
  }
  /** Draw a random number from the discrete distribution.
   * \return Sampled value
   */
  int operator()() { return distribution(engine); }

 private:
  /** The distribution object that is being used. */
  std::discrete_distribution<> distribution;
};

/**
 * Draws a random number from a Cauchy distribution (sometimes also called
 * Lorentz or non-relativistic Breit-Wigner distribution) with the given
 * parameters (constant width!) inside the range [min,max]. This function is
 * similar to std::cauchy_distribution, but can return values inside a limited
 * interval.
 * \param pole Pole parameter of the Cauchy function, i.e. location of the peak.
 * \param width Width parameter of the Cauchy function, determining the
 * sharpness of the peak.
 * \param min Minimum value to be returned.
 * \param max Maximum value to be returned.
 * \return Sampled random number.
 */
template <typename T = double>
T cauchy(T pole, T width, T min, T max) {
  /* Use double-precision variables, in order to work around a glibc bug in
   * tanf:
   * https://sourceware.org/bugzilla/show_bug.cgi?id=18221 */
  const double u_min = std::atan((min - pole) / width);
  const double u_max = std::atan((max - pole) / width);
  const double u = uniform(u_min, u_max);
  return pole + width * std::tan(u);
}

/**
 * Draws a random number from a beta-distribution, where probability density of
 * \f$x\f$ is \f$p(x) = frac{\Gamma(a)\Gamma(b)}{Gamma(a+b)}
 *  x^{a-1} (1-x)^{b-1}\f$. This distribution is necessary for string
 *  formation. The implementation uses a property connecting beta distribution
 * to gamma-distribution. Interchanging a and b will not change results.
 *
 * \param a Shape parameter.
 * \param b Scale parameter.
 * \return Sampled random number.
 */
template <typename T = double>
T beta(T a, T b) {
  // Otherwise the integral over probability density diverges
  assert(a > T(0.0) && b > T(0.0));
  const T x1 = std::gamma_distribution<T>(a)(engine);
  const T x2 = std::gamma_distribution<T>(b)(engine);
  return x1 / (x1 + x2);
}

/**
 * Draws a random number from a beta-distribution with a = 0. In this case
 * the probability density is \f$p(x) = 1/x (1-x)^b\f$. The integral from
 * 0 to 1 over this distribution diverges, so the sampling is performed
 * in the interval (xmin, 1). This distribution is necessary for string
 * formation. The implementation uses the following property:
 * \f$p(x)dx = dx/x (1-x)^b = (1-x)^b d ln(x) = (1 - e^{-y})^b dy\f$, where
 * \f$ y = - ln(x) \f$.
 *
 * \param xmin Minimal sampled value.
 * \param b Second shape parameter.
 * \return Sampled random number.
 */
template <typename T = double>
T beta_a0(T xmin, T b) {
  assert(xmin > T(0.0) && xmin < T(1.0));
  T y;
  do {
    y = uniform(0.0, -std::log(xmin));
  } while (std::pow((1.0 - std::exp(-y)), b) < canonical());
  return std::exp(-y);
}

}  // namespace random
}  // namespace smash

#endif  // SRC_INCLUDE_RANDOM_H_
