/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_INTEGRATE_H_
#define SRC_INCLUDE_INTEGRATE_H_

#include <cuba.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>

#include "cxx14compat.h"
#include "fpenvironment.h"
#include "random.h"

namespace smash {

/**
 * A deleter type for std::unique_ptr to be used with
 * `gsl_integration_workspace` pointers. This will call
 * `gsl_integration_workspace_free` instead of `delete`.
 */
struct GslWorkspaceDeleter {
  /// The class has no members, so this is a noop.
  constexpr GslWorkspaceDeleter() = default;

  /// Frees the gsl_integration_cquad_workspace resource if it is non-zero.
  void operator()(gsl_integration_cquad_workspace *ptr) const {
    if (ptr == nullptr) {
      return;
    }
    gsl_integration_cquad_workspace_free(ptr);
  }
};

/**
 * The result type returned from integrations,
 * containing the value and an error.
 */
class Result : public std::pair<double, double> {
  /// The data type to store the value and the error of the integration
  using Base = std::pair<double, double>;

 public:
  /// Forward the pair constructors.
  using Base::pair;

  /// Conversion to double yields the value of the integral.
  operator double() const { return Base::first; }

  /// Access the first entry in the pair as the value.
  double value() const { return Base::first; }

  /// Access the second entry in the pair as the absolute error.
  double error() const { return Base::second; }

  /**
   * Check whether the error is small and alert if it is not.
   *
   * \param[in] integration_name Name of the integration used
   *            for error message.
   * \param[in] relative_tolerance Relative error tolerance.
   * \param[in] absolute_tolerance Absolute error tolerance.
   */
  void check_error(const std::string &integration_name,
                   double relative_tolerance = 5e-4,
                   double absolute_tolerance = 1e-9) const {
    const double allowed_error =
        std::max(absolute_tolerance, value() * relative_tolerance);
    if (error() > allowed_error) {
      std::stringstream error_msg;
      error_msg << integration_name << " resulted in I = " << value() << " ± "
                << error()
                << ", but the required precision is either absolute error < "
                << absolute_tolerance << " or relative error < "
                << relative_tolerance << std::endl;
      throw std::runtime_error(error_msg.str());
    }
  }
};

/**
 * A C++ interface for numerical integration in one dimension
 * with the GSL CQUAD integration functions.
 *
 * Example:
 * \code
 * Integrator integrate;
 * const auto result = integrate(0.1, 0.9, [](double x) { return x * x; });
 * \endcode
 */
class Integrator {
 public:
  /**
   * Construct an integration functor with the given \p workspace_size.
   *
   * \note Since the workspace is allocated in the constructor and deallocated
   * on destruction, you should not recreate Integrator objects unless required.
   * Thus, if you want to calculate multiple integrals with the same \p
   * workspace_size, keep the Integrator object around.
   *
   * \param workspace_size The internal workspace is allocated such that it can
   *                       hold the given number of double precision intervals,
   *                       their integration results, and error estimates.
   *                       It also determines the maximum number of subintervals
   *                       the integration algorithm will use.
   */
  explicit Integrator(size_t workspace_size = 1000)
      : workspace_(gsl_integration_cquad_workspace_alloc(workspace_size)) {}

  /**
   * The function call operator implements the integration functionality.
   *
   * \param[in] a The lower limit of the integral.
   * \param[in] b The upper limit of the integral.
   * \tparam F Type of the integrand function.
   * \param[in] fun The callable to integrate over. This callable may be a
   *            function pointer, lambda, or a functor object. In any case, the
   *            callable must return a `double` and take a single `double`
   *            argument. If you want to pass additional data to the callable
   *            you can e.g. use lambda captures.
   * \return Pair of integral value and absolute error estimate.
   */
  template <typename F>
  Result operator()(double a, double b, F &&fun) {
    Result result;
    const gsl_function gslfun{
        /* important! The lambda cannot use captures, otherwise the
         * conversion to a function pointer type is impossible. */
        [](double x, void *type_erased) -> double {
          auto &&f = *static_cast<F *>(type_erased);
          return f(x);
        },
        &fun};
    // We disable float traps when calling GSL code we cannot control.
    DisableFloatTraps guard;
    const int error_code = gsl_integration_cquad(
        &gslfun, a, b, accuracy_absolute_, accuracy_relative_, workspace_.get(),
        &result.first, &result.second,
        nullptr /* Don't store the number of evaluations */);
    if (error_code) {
      std::stringstream err;
      err << "GSL 1D deterministic integration: " << gsl_strerror(error_code);
      throw std::runtime_error(err.str());
    }
    result.check_error("GSL 1D deterministic integration", accuracy_relative_,
                       accuracy_absolute_);
    return result;
  }

 private:
  /// Holds the workspace pointer.
  std::unique_ptr<gsl_integration_cquad_workspace, GslWorkspaceDeleter>
      workspace_;

  /// Parameter to the GSL integration function: desired absolute error limit
  const double accuracy_absolute_ = 1.0e-9;

  /// Parameter to the GSL integration function: desired relative error limit
  const double accuracy_relative_ = 5.0e-4;
};

/**
 * A C++ interface for numerical integration in one dimension
 * with the GSL Monte-Carlo integration functions.
 *
 * Example:
 * \code
 * Integrator integrate;
 * const auto result = integrate(0.1, 0.9,
 *                               [](double x) { return x * x; });
 * \endcode
 */
class Integrator1dMonte {
 public:
  /**
   * Construct an integration functor.
   *
   * \param[in] num_calls The desired number of calls to the integrand function
   *                  (defaults to 1E6 if omitted), i.e. how often the integrand
   *                  is sampled in the Monte-Carlo integration. Larger numbers
   *                  lead to a more precise result, but also to increased
   *                  runtime.
   *
   * \note Since the workspace is allocated in the constructor and deallocated
   * on destruction, you should not recreate Integrator objects unless required.
   * Thus, if you want to calculate multiple integrals with the same \p
   * workspace_size, keep the Integrator object around.
   */
  explicit Integrator1dMonte(size_t num_calls = 1E6)
      : state_(gsl_monte_plain_alloc(1)),
        rng_(gsl_rng_alloc(gsl_rng_mt19937)),
        number_of_calls_(num_calls) {
    gsl_monte_plain_init(state_);
    // initialize the GSL RNG with a random seed
    const uint32_t seed = random::uniform_int(0ul, ULONG_MAX);
    gsl_rng_set(rng_, seed);
  }

  /// Destructor: Clean up internal state and RNG.
  ~Integrator1dMonte() {
    gsl_monte_plain_free(state_);
    gsl_rng_free(rng_);
  }

  /**
   * The function call operator implements the integration functionality.
   *
   * \param[in] min The lower limit of the integration.
   * \param[in] max The upper limit of the integration.
   * \tparam F Type of the integrand function.
   * \param[in] fun
   *            The callable to integrate over. This callable may be a function
   *            pointer, lambda, or a functor object. In any case, the callable
   *            must return a `double` and take two `double` arguments. If you
   *            want to pass additional data to the callable you can e.g. use
   *            lambda captures.
   * \return Pair of integral value and absolute error estimate.
   */
  template <typename F>
  Result operator()(double min, double max, F &&fun) {
    Result result = {0, 0};

    const double lower[1] = {min};
    const double upper[1] = {max};

    if (max <= min)
      return result;

    const gsl_monte_function monte_fun{
        // trick: pass integrand function as 'params'
        [](double *x, size_t /*dim*/, void *params) -> double {
          auto &&f = *static_cast<F *>(params);
          return f(x[0]);
        },
        1, &fun};

    const int error_code =
        gsl_monte_plain_integrate(&monte_fun, lower, upper, 1, number_of_calls_,
                                  rng_, state_, &result.first, &result.second);
    if (error_code) {
      std::stringstream err;
      err << "GSL 1D Monte-Carlo integration: " << gsl_strerror(error_code);
      throw std::runtime_error(err.str());
    }

    result.check_error("GSL 1D Monte-Carlo integration");

    return result;
  }

 private:
  /// internal state of the Monte-Carlo integrator
  gsl_monte_plain_state *state_;

  /// random number generator
  gsl_rng *rng_;

  /// number of calls to the integrand
  const std::size_t number_of_calls_;
};

/**
 * A C++ interface for numerical integration in two dimensions
 * with the GSL Monte-Carlo integration functions.
 *
 * Example:
 * \code
 * Integrator integrate;
 * const auto result = integrate(0.1, 0.9, 0., 0.5,
 *                               [](double x, double y) { return x * y; });
 * \endcode
 */
class Integrator2d {
 public:
  /**
   * Construct an integration functor.
   *
   * \param[in] num_calls The desired number of calls to the integrand function
   *                  (defaults to 1E6 if omitted), i.e. how often the integrand
   *                  is sampled in the Monte-Carlo integration. Larger numbers
   *                  lead to a more precise result, but also to increased
   *                  runtime.
   *
   * \note Since the workspace is allocated in the constructor and deallocated
   * on destruction, you should not recreate Integrator objects unless required.
   * Thus, if you want to calculate multiple integrals with the same \p
   * workspace_size, keep the Integrator object around.
   */
  explicit Integrator2d(size_t num_calls = 1E6)
      : state_(gsl_monte_plain_alloc(2)),
        rng_(gsl_rng_alloc(gsl_rng_mt19937)),
        number_of_calls_(num_calls) {
    gsl_monte_plain_init(state_);
    // initialize the GSL RNG with a random seed
    const uint32_t seed = random::uniform_int(0ul, ULONG_MAX);
    gsl_rng_set(rng_, seed);
  }

  /// Destructor: Clean up internal state and RNG.
  ~Integrator2d() {
    gsl_monte_plain_free(state_);
    gsl_rng_free(rng_);
  }

  /**
   * The function call operator implements the integration functionality.
   *
   * \param[in] min1 The lower limit in the first dimension.
   * \param[in] max1 The upper limit in the first dimension.
   * \param[in] min2 The lower limit in the second dimension.
   * \param[in] max2 The upper limit in the second dimension.
   * \tparam F Type of the integrand function.
   * \param[in] fun
   *            The callable to integrate over. This callable may be a function
   *            pointer, lambda, or a functor object. In any case, the callable
   *            must return a `double` and take two `double` arguments. If you
   *            want to pass additional data to the callable you can e.g. use
   *            lambda captures.
   * \return Pair of integral value and absolute error estimate.
   */
  template <typename F>
  Result operator()(double min1, double max1, double min2, double max2,
                    F &&fun) {
    Result result = {0, 0};

    const double lower[2] = {min1, min2};
    const double upper[2] = {max1, max2};

    if (max1 <= min1 || max2 <= min2)
      return result;

    const gsl_monte_function monte_fun{
        // trick: pass integrand function as 'params'
        [](double *x, size_t /*dim*/, void *params) -> double {
          auto &&f = *static_cast<F *>(params);
          return f(x[0], x[1]);
        },
        2, &fun};

    gsl_monte_plain_integrate(&monte_fun, lower, upper, 2, number_of_calls_,
                              rng_, state_, &result.first, &result.second);

    return result;
  }

 private:
  /// internal state of the Monte-Carlo integrator
  gsl_monte_plain_state *state_;

  /// random number generator
  gsl_rng *rng_;

  /// number of calls to the integrand
  const std::size_t number_of_calls_;
};

/**
 * This is a wrapper for the integrand, so we can pass the limits as well for
 * renormalizing to the unit cube.
 *
 * \tparam F Type of the integrand function.
 */
template <typename F>
struct Integrand2d {
  /// the lower bound of the first integrated variable
  double min1;
  /// the integration range of the first integrated variable
  double diff1;
  /// the lower bound of the second integrated variable
  double min2;
  /// the integration range of the second integrated variable
  double diff2;
  /// the integrated function
  F f;
};

/**
 * A C++ interface for numerical integration in two dimensions with the Cuba
 * Cuhre integration function.
 *
 * The algorithm is deterministic and well-suited for low dimensions, where it
 * can reach good accuracy.
 *
 * Example:
 * \code
 * Integrator2dCuhre integrate;
 * const auto result = integrate(0.1, 0.9, 0., 0.5,
 *                               [](double x, double y) { return x * y; });
 * \endcode
 */
class Integrator2dCuhre {
 public:
  /**
   * Construct an integration functor.
   *
   * \param[in] num_calls The maximum number of calls to the integrand function
   *                  (defaults to 1E6 if omitted), i.e. how often the integrand
   *                  can be sampled in the integration. Larger numbers lead to
   * a
   *                  more precise result, but possibly to increased runtime.
   * \param[in] epsrel    The desired relative accuracy (1E-3 by default).
   * \param[in] epsabs    The desired absolute accuracy (1E-3 by default).
   */
  explicit Integrator2dCuhre(int num_calls = 1e6, double epsrel = 5e-4,
                             double epsabs = 1e-9)
      : maxeval_(num_calls), epsrel_(epsrel), epsabs_(epsabs) {}

  /**
   * The function call operator implements the integration functionality.
   *
   * \param[in] min1 The lower limit in the first dimension.
   * \param[in] max1 The upper limit in the first dimension.
   * \param[in] min2 The lower limit in the second dimension.
   * \param[in] max2 The upper limit in the second dimension.
   * \tparam F Type of the integrand function.
   * \param[in] fun
   *            The callable to integrate over. This callable may be a function
   *            pointer, lambda, or a functor object. In any case, the callable
   *            must return a `double` and take two `double` arguments. If you
   *            want to pass additional data to the callable you can e.g. use
   *            lambda captures.
   * \return Pair of integral value and absolute error estimate.
   */
  template <typename F>
  Result operator()(double min1, double max1, double min2, double max2, F fun) {
    Result result = {0., 0.};

    if (max1 < min1 || max2 < min2) {
      bool tolerable = (max1 - min1 > -1.e-16) && (max2 - min2 > -1.e-16);
      if (tolerable) {
        return result;
      }
      std::stringstream err;
      err << "Integrator2dCuhre got wrong integration limits: [" << min1 << ", "
          << max1 << "], [" << min2 << ", " << max2 << "]";
      throw std::invalid_argument(err.str());
    }

    Integrand2d<F> f_with_limits = {min1, max1 - min1, min2, max2 - min2, fun};

    const integrand_t cuhre_fun{[](const int * /* ndim */, const cubareal xx[],
                                   const int * /* ncomp */, cubareal ff[],
                                   void *userdata) -> int {
      auto i = static_cast<Integrand2d<F> *>(userdata);
      /* We have to transform the integrand to the unit cube.
       * This is what Cuba expects. */
      ff[0] = (i->f)(i->min1 + i->diff1 * xx[0], i->min2 + i->diff2 * xx[1]) *
              i->diff1 * i->diff2;
      return 0;
    }};

    const int ndim = 2;
    const int ncomp = 1;
    void *userdata = &f_with_limits;
    const int nvec = 1;
    const int flags = 0;  // Use the defaults.
    const int mineval = 0;
    const int maxeval = maxeval_;
    const int key = -1;  // Use the default.
    const char *statefile = nullptr;
    void *spin = nullptr;

    Cuhre(ndim, ncomp, cuhre_fun, userdata, nvec, epsrel_, epsabs_, flags,
          mineval, maxeval, key, statefile, spin, &nregions_, &neval_, &fail_,
          &result.first, &result.second, &prob_);

    if (fail_) {
      std::stringstream err;
      err << "After " << neval_ << " evaluations "
          << "Cuhre integration from Cuba reports error code " << fail_;
      throw std::runtime_error(err.str());
    }
    result.check_error("Cuba integration ", epsrel_, epsabs_);

    return result;
  }

 private:
  /// The (approximate) maximum number of integrand evaluations allowed.
  int maxeval_;
  /// Requested relative accuracy.
  double epsrel_;
  /// Requested absolute accuracy.
  double epsabs_;
  /// Actual number of subregions needed.
  int nregions_;
  /// Actual number of integrand evaluations needed.
  int neval_;
  /** 
   * An error flag.
   *
   * 0 if the desired accuracy was reached, -1 if the dimension is out of
   * range, larger than 0 if the accuracy goal was not met within the maximum
   * number of evaluations.
   */
  int fail_;
  /**
   * The chi^2 probability that the error is not a reliable estimate of the
   * true integration error.
   */
  double prob_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_INTEGRATE_H_
