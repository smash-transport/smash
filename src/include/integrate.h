/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_INTEGRATE_H_
#define SRC_INCLUDE_INTEGRATE_H_

#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_plain.h>
#include <tuple>

#include "cxx14compat.h"

namespace Smash {

/**
 * A deleter type for std::unique_ptr to be used with
 * `gsl_integration_workspace` pointers. This will call
 * `gsl_integration_workspace_free` instead of `delete`.
 */
struct GslWorkspaceDeleter {
  /// The class has no members, so this is a noop.
  constexpr GslWorkspaceDeleter() = default;

  /// frees the gsl_integration_workspace resource if it is non-zero.
  void operator()(gsl_integration_workspace *ptr) const {
    if (ptr == nullptr) {
      return;
    }
    gsl_integration_workspace_free(ptr);
  }
};


/** The result type returned from integrations,
 * containing the value and an error. */
class Result : public std::pair<double, double> {
  using Base = std::pair<double, double>;

 public:
  /// forward the pair constructors
  using Base::pair;

  /// conversion to double yields the value of the integral
  operator double() const { return Base::first; }

  /// access the first entry in the pair as the value
  double value() const { return Base::first; }

  /// access the second entry in the pair as the absolute error
  double error() const { return Base::second; }
};


/**
 * A C++ interface for numerical integration in one dimension
 * with the GSL integration functions.
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
   */
  explicit Integrator(int workspace_size)
      : workspace_(gsl_integration_workspace_alloc(workspace_size)) {}

  /// Convenience overload of the above with a workspace size of 1000.
  Integrator() : Integrator(1000) {}

  /**
   * The function call operator implements the integration functionality.
   *
   * \param a The lower limit of the integral.
   * \param b The upper limit of the integral.
   * \param fun The callable to integrate over. This callable may be a function
   *            pointer, lambda, or a functor object. In any case, the callable
   *            must return a `double` and take a single `double` argument. If
   *            you want to pass additional data to the callable you can e.g.
   *            use lambda captures.
   */
  template <typename F>
  Result operator()(double a, double b, F &&fun) {
    Result result;
    const gsl_function gslfun{
        // important! The lambda cannot use captures, otherwise the
        // conversion to a function pointer type is impossible.
        [](double x, void *type_erased) -> double {
          auto &&f = *static_cast<F *>(type_erased);
          return f(x);
        },
        &fun};
    gsl_integration_qag(&gslfun, a, b,
                        accuracy_absolute_,  // epsabs
                        accuracy_relative_,  // epsrel
                        subintervals_max_,   // limit
                        gauss_points_,       // key
                        workspace_.get(), &result.first, &result.second);
    return result;
  }

 private:
  /// Holds the workspace pointer.
  std::unique_ptr<gsl_integration_workspace, GslWorkspaceDeleter> workspace_;

  /// Parameter to the GSL integration function: desired absolute error limit
  double accuracy_absolute_ = 1.0e-5;

  /// Parameter to the GSL integration function: desired relative error limit
  double accuracy_relative_ = 5.0e-4;

  /// Parameter to the GSL integration function: maximum number of subintervals
  /// (may not exceed workspace size)
  std::size_t subintervals_max_ = 500;

  /// Parameter to the GSL integration function: integration rule
  int gauss_points_ = GSL_INTEG_GAUSS21;
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
   * \param num_calls The desired number of calls to the integrand function
   *                  (defaults to 1E6 if omitted), i.e. how often the integrand
   *                  is sampled in the Monte-Carlo integration. Larger numbers
   *                  lead to a more precise result, but also to increased runtime.
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
  }

  /**
   * Destructor: Clean up internal state and RNG.
   */
  ~Integrator2d() {
    gsl_monte_plain_free(state_);
    gsl_rng_free(rng_);
  }

  /**
   * The function call operator implements the integration functionality.
   *
   * \param min1 The lower limit in the first dimension.
   * \param max1 The upper limit in the first dimension.
   * \param min2 The lower limit in the second dimension.
   * \param max2 The upper limit in the second dimension.
   * \param fun The callable to integrate over. This callable may be a function
   *            pointer, lambda, or a functor object. In any case, the callable
   *            must return a `double` and take two `double` arguments. If you
   *            want to pass additional data to the callable you can e.g. use
   *            lambda captures.
   */
  template <typename F>
  Result operator()(double min1, double max1,
                    double min2, double max2, F &&fun) {
    Result result = {0, 0};

    const double lower[2] = {min1, min2};
    const double upper[2] = {max1, max2};

    if (max1 <= min1 || max2 <= min2)
      return result;

    const gsl_monte_function monte_fun {
        // trick: pass integrand function as 'params'
        [](double *x, size_t dim, void *params) -> double {
          auto &&f = *static_cast<F *>(params);
          assert(dim == 2);
          return f(x[0], x[1]);
        },
        2, &fun
    };

    gsl_monte_plain_integrate(&monte_fun, lower, upper, 2,
                              number_of_calls_, rng_, state_,
                              &result.first, &result.second);

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

}  // namespace Smash

#endif  // SRC_INCLUDE_INTEGRATE_H_
