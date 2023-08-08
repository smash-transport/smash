/*
 *
 *    Copyright (c) 2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_ROOTSOLVER_H_
#define SRC_INCLUDE_SMASH_ROOTSOLVER_H_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <functional>
#include <memory>
#include <string>

#include "logging.h"

namespace smash {
static constexpr int LRootSolver = LogArea::RootSolver::id;
/**
 * A class used for calculating the root of a one-dimensional equation.
 * It takes care of all tecnicalities and provides an interace to GSL.
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
   * The function of which a root should be found in the form that GSL expects
   *
   * \param[in] x Position on x-axis where the function should be evaluated
   *
   * \return Feedback whether the equation was evaluatd successfully
   */
  static double gsl_func(const double x, void *) { return (*root_eq_)(x); }

  /**
   * Attempt to find a root in a given interval
   *
   * \param[in] initial_guess_low Lower boundary of the interval
   * \param[in] initial_guess_high Higher boundary of the interval
   * \param[in] itermax maximum number of steps for root finding
   * \param[out] root Root of the function if found successfully
   * \return true if a root was successfully found
   */
  bool try_find_root(double initial_guess_low, double initial_guess_high,
                     size_t itermax, double &root) {
    // check if root is in the given interval
    if ((*root_eq_)(initial_guess_low) * (*root_eq_)(initial_guess_high) > 0) {
      return false;
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
        root = 0.5 * (xlow + xhigh);
        gsl_root_fsolver_free(Root_finder_);
        return true;
      }
    } while (status == GSL_CONTINUE && iter < itermax);
    gsl_root_fsolver_free(Root_finder_);
    return false;
  }

 private:
  /// GSL solver to use for rootfingding
  const gsl_root_fsolver_type *Solver_name_ = gsl_root_fsolver_brent;

  /// GSL rootfinding object to take care of root findung
  gsl_root_fsolver *Root_finder_;

  /// Static pointer to the function to solve
  static inline std::unique_ptr<std::function<double(double)>> root_eq_ =
      nullptr;

  /// Expected precision of the root
  double solution_precision_ = 1e-7;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_ROOTSOLVER_H_
