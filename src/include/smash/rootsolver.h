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

#include <functional>
#include <memory>
#include <string>

#include "gsl/gsl_multiroots.h"
#include "gsl/gsl_vector.h"

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
   * \param[inout] roots_array C-Array with x-values of the function
   * \param[inout] function c-Array with f(x)-values of the function
   * \return Feedback whether the equation was evaluatd successfully
   */
  static int gsl_func(const gsl_vector *roots_array, void *,
                      gsl_vector *function) {
    double x = gsl_vector_get(roots_array, 0);
    gsl_vector_set(function, 0, (*root_eq_)(x));
    return GSL_SUCCESS;
  }

  /**
   * Attempt to find a root with a given initial guess
   *
   * \param[in] initial_guess First x-value to start root finding
   * \param[in] itermax maximum number of steps for root finding
   * \param[out] root
   * \return true if a root was successfully found
   */
  bool try_find_root(double initial_guess, size_t itermax, double &root) {
    gsl_multiroot_function function_GSL = {&(gsl_func), 1, nullptr};
    int status = GSL_CONTINUE;
    size_t iter = 0;
    gsl_vector *roots_array = gsl_vector_alloc(1);
    Root_finder_ = gsl_multiroot_fsolver_alloc(Solver_name_, 1);
    gsl_vector_set(roots_array, 0, initial_guess);
    gsl_multiroot_fsolver_set(Root_finder_, &function_GSL, roots_array);
    do {
      iter++;
      status = gsl_multiroot_fsolver_iterate(Root_finder_);
      if (status != GSL_SUCCESS) {
        logg[LRootSolver].debug("GSL ERROR in root finding: " +
                                static_cast<std::string>(gsl_strerror(status)) +
                                "\n with starting value " +
                                std::to_string(initial_guess));
        break;
      }
      status =
          gsl_multiroot_test_residual(Root_finder_->f, solution_precision_);
      if (status == GSL_SUCCESS) {
        root = gsl_vector_get(Root_finder_->x, 0);
        gsl_multiroot_fsolver_free(Root_finder_);
        gsl_vector_free(roots_array);
        return true;
      }
    } while (status == GSL_CONTINUE && iter < itermax);
    gsl_multiroot_fsolver_free(Root_finder_);
    gsl_vector_free(roots_array);
    return false;
  }

 private:
  /// GSL solver to use for rootfingding
  const gsl_multiroot_fsolver_type *Solver_name_ =
      gsl_multiroot_fsolver_hybrids;

  /// GSL rootfinding object to take care of root findung
  gsl_multiroot_fsolver *Root_finder_;

  /// Static pointer to the function to solve
  static inline std::unique_ptr<std::function<double(double)>> root_eq_ =
      nullptr;

  /// Expected precision of the root
  double solution_precision_ = 1e-7;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_ROOTSOLVER_H_

