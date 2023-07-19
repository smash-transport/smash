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

#include "gsl/gsl_multiroots.h"
#include "gsl/gsl_vector.h"

#include "logging.h"

namespace smash {
static constexpr int LRootSolver = LogArea::RootSolver::id;
class RootSolver1D {
 public:
  RootSolver1D(std::function<double(double)> eq) {
    root_eq_ = std::make_unique<std::function<double(double)>>(eq);
  }

  ~RootSolver1D() { root_eq_ = nullptr; }

  static int gsl_func(const gsl_vector *roots_array, void *,
                      gsl_vector *function) {
    double x = gsl_vector_get(roots_array, 0);
    gsl_vector_set(function, 0, (*root_eq_)(x));
    return GSL_SUCCESS;
  }

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
  const gsl_multiroot_fsolver_type *Solver_name_ =
      gsl_multiroot_fsolver_hybrids;

  gsl_multiroot_fsolver *Root_finder_;

  static inline std::unique_ptr<std::function<double(double)>> root_eq_ =
      nullptr;

  double solution_precision_ = 1e-7;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_ROOTSOLVER_H_