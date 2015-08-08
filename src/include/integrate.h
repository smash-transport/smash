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

/**
 * A C++ interface to GSL integration functions.
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

  /// The result type returned from integrations containing the value and an
  /// error.
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

}  // namespace Smash

#endif  // SRC_INCLUDE_INTEGRATE_H_
