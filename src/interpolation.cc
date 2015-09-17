/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/interpolation.h"

#include <iostream>
#include <sstream>
#include <tuple>

static std::pair<std::vector<double>, std::vector<double>> sort_data(
        const std::vector<double>& x, const std::vector<double> y) {
    const auto p = generate_sort_permutation(
        x, [&](double const& a, double const& b) {
          return a < b;
        });
    const std::vector<double> sorted_x = std::move(apply_permutation(x, p));
    const std::vector<double> sorted_y = std::move(apply_permutation(y, p));
    for (size_t i = 0; i < sorted_x.size() - 1; i++) {
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wfloat-equal"
        if (sorted_x[i] == sorted_x[i + 1]) {
        #pragma GCC diagnostic pop
          std::stringstream error_msg;
          error_msg << "InterpolateDataSpline: Each x value must be unique. \""
                    << sorted_x[i] << "\" was found twice.";
          throw std::runtime_error(error_msg.str());
        }
    }
    return std::make_pair(std::move(sorted_x), std::move(sorted_y));
}

InterpolateDataSpline::InterpolateDataSpline(const std::vector<double>& x,
                                             const std::vector<double>& y) {
    const auto N = x.size();
    if (y.size() != N) {
      throw std::runtime_error("Need two vectors of equal length "
                               "for interpolation.");
    }
    if (N < 3) {
      throw std::runtime_error("Need at least 3 data points "
                               "for cubic spline interpolation.");
    }
    std::vector<double> sorted_x, sorted_y;
    std::tie(sorted_x, sorted_y) = sort_data(x, y);
    first_x_ = sorted_x.front();
    last_x_ = sorted_x.back();
    first_y_ = sorted_y.front();
    last_y_ = sorted_y.back();
    acc_ = gsl_interp_accel_alloc();
    spline_ = gsl_spline_alloc(gsl_interp_cspline, N);
    gsl_spline_init(spline_, &(*sorted_x.begin()), &(*sorted_y.begin()), N);
}

InterpolateDataSpline::~InterpolateDataSpline() {
    gsl_spline_free(spline_);
    gsl_interp_accel_free(acc_);
}

double InterpolateDataSpline::operator()(double xi) const {
    // constant extrapolation
    if (xi < first_x_) {
        return first_y_;
    }
    if (xi > last_x_) {
        return last_y_;
    }
    // cubic spline interpolation
    return gsl_spline_eval(spline_, xi, acc_);
}

InterpolateDataBSpline::InterpolateDataBSpline(const std::vector<double>& x,
                                               const std::vector<double>& y,
                                               size_t nbreak) {
    const size_t N = x.size();
    if (y.size() != N) {
      throw std::runtime_error("Need two vectors of equal length "
                               "for interpolation.");
    }
    // Make sure the data is sorted and unambigous.
    std::vector<double> sorted_x, sorted_y;
    std::tie(sorted_x, sorted_y) = sort_data(x, y);

    // Allocate the B-spline and fit workspace.
    constexpr size_t k = 4;
    spline_workspace_ = gsl_bspline_alloc(k, nbreak);
    const int err = gsl_bspline_knots_uniform(sorted_x.front(), sorted_x.back(),
                                              spline_workspace_);
    if (err) {
      throw std::runtime_error("could not construct B-spline");
    }
    const size_t ncoeffs = gsl_bspline_ncoeffs(spline_workspace_);
    B_ = gsl_vector_alloc(ncoeffs);
    auto X = gsl_matrix_alloc(N, ncoeffs);
    c_ = gsl_vector_alloc(ncoeffs);
    cov_ = gsl_matrix_alloc(ncoeffs, ncoeffs);
    fit_workspace_ = gsl_multifit_linear_alloc(N, ncoeffs);

    // Construct the fit matrix X.
    for (size_t i = 0; i < N; i++) {
        const double xi = sorted_x[i];
        // Compute B_j(xi) for all j.
        gsl_bspline_eval(xi, B_, spline_workspace_);
        // Fill in row i of X.
        for (size_t j = 0; j < ncoeffs; j++) {
            const double Bj = gsl_vector_get(B_, j);
            gsl_matrix_set(X, i, j, Bj);
        }
    }

    // Do the fit.
    double chisq;
    gsl_vector gsl_y;
    gsl_y.size = sorted_y.size();
    gsl_y.stride = 1;
    gsl_y.data = &sorted_y.front();
    gsl_y.block = nullptr;
    gsl_y.owner = 0;
    gsl_multifit_linear(X, &gsl_y, c_, cov_, &chisq, fit_workspace_);

    gsl_matrix_free(X);
}

InterpolateDataBSpline::~InterpolateDataBSpline() {
    gsl_bspline_free(spline_workspace_);
    gsl_multifit_linear_free(fit_workspace_);
    gsl_vector_free(B_);
    gsl_vector_free(c_);
    gsl_matrix_free(cov_);
}

double InterpolateDataBSpline::operator()(double xi) const {
    double yi;
    double yerr;
    gsl_bspline_eval(xi, B_, spline_workspace_);
    gsl_multifit_linear_est(B_, c_, cov_, &yi, &yerr);
    return yi;
}
