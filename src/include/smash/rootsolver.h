/*
 *
 *    Copyright (c) 2023,2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_ROOTSOLVER_H_
#define SRC_INCLUDE_SMASH_ROOTSOLVER_H_

#include <functional>
#include <memory>
#include <optional>
#include <string>

#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_roots.h"

#include "logging.h"

namespace smash {
static constexpr int LRootSolver = LogArea::RootSolver::id;
/**
 * A class used for calculating the root of a one-dimensional equation.
 * It takes care of all technicalities and provides an interface to GSL.
 */
class RootSolver1D {
 public:
  /**
   * Construct a new Root Solver 1D object
   *
   * \param[in] eq The function of which a root is desired
   */
  explicit RootSolver1D(std::function<double(double)> eq) {
    root_eq_ = std::make_unique<std::function<double(double)>>(eq);
  }

  ~RootSolver1D() { root_eq_ = nullptr; }

  /**
   * Attempt to find a root in a given interval
   *
   * \param[in] initial_guess_low Lower boundary of the interval
   * \param[in] initial_guess_high Higher boundary of the interval
   * \param[in] itermax maximum number of steps for root finding
   * \return the root if it was found
   */
  std::optional<double> try_find_root(double initial_guess_low,
                                      double initial_guess_high,
                                      size_t itermax) {
    // check if root is in the given interval
    if ((*root_eq_)(initial_guess_low) * (*root_eq_)(initial_guess_high) > 0) {
      logg[LRootSolver].trace()
          << "Function has same sign at both ends of the interval ["
          << initial_guess_low << ", " << initial_guess_high
          << "]. Root can't be found in this interval.";
      return std::nullopt;
    }
    gsl_function function_GSL = {&(gsl_func), nullptr};
    int status = GSL_CONTINUE;
    size_t iter = 0;
    Root_finder_ = gsl_root_fsolver_alloc(Solver_name_);
    gsl_root_fsolver_set(Root_finder_, &function_GSL, initial_guess_low,
                         initial_guess_high);
    do {
      iter++;
      status = gsl_root_fsolver_iterate(Root_finder_);
      if (status != GSL_SUCCESS) {
        logg[LRootSolver].debug("GSL ERROR in root finding: " +
                                static_cast<std::string>(gsl_strerror(status)));
        break;
      }
      double xlow = gsl_root_fsolver_x_lower(Root_finder_);
      double xhigh = gsl_root_fsolver_x_upper(Root_finder_);
      status = gsl_root_test_interval(xlow, xhigh, 0, solution_precision_);
      if (status == GSL_SUCCESS) {
        double root = 0.5 * (xlow + xhigh);
        gsl_root_fsolver_free(Root_finder_);
        return root;
      }
    } while (status == GSL_CONTINUE && iter < itermax);
    gsl_root_fsolver_free(Root_finder_);
    return std::nullopt;
  }

 private:
  /// GSL solver to use for root finding
  const gsl_root_fsolver_type *Solver_name_ = gsl_root_fsolver_brent;

  /// GSL root finding object to take care of root finding
  gsl_root_fsolver *Root_finder_ = nullptr;

  /** Static pointer to the function to solve
   * \note This member has to be \c static since it is used inside the `static
   * gsl_func` method. As \c static methods exist and can be used without a
   * class instance, it is not possible to use non-static members inside them
   * (non-static members \b are bound to a class instance). \see gsl_func to
   * understand why that method needs to be \c static as well.
   */
  static inline std::unique_ptr<std::function<double(double)>> root_eq_ =
      nullptr;

  /// Expected precision of the root
  double solution_precision_ = 1e-7;

  /**
   * The function of which a root should be found in the form that GSL expects
   *
   * \param[in] x Argument of the function
   *
   * \return The value of the function for the given argument x
   *
   * \note
   * This function needs to be \c static because GSL uses an object of type \c
   * gsl_function that has to be initialised with a pointer to a function and
   * this cannot be created from a class non-static method as a class non-static
   * method cannot exist as entity without an instance of that type. On the
   * contrary, a pointer to a static member function can be used as a normal
   * pointer to a free function, because a \c static method exists without
   * needing any instance of the class.
   */
  static double gsl_func(const double x, void *) { return (*root_eq_)(x); }
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_ROOTSOLVER_H_
