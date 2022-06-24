/*
 *
 *    Copyright (c) 2020,2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_INTERPOLATION2D_H_
#define SRC_INCLUDE_SMASH_INTERPOLATION2D_H_

#include <gsl/gsl_spline2d.h>

#include <vector>

namespace smash {

/// Represent a bicubic spline interpolation.
class InterpolateData2DSpline {
 public:
  /**
   * Interpolate function f given discrete samples f(x_i, y_i) = z_i.
   *
   * \param x x-values.
   * \param y y-values.
   * \param z z-values
   * \return The interpolation function.
   *
   * A bicubic spline interpolation is used.
   * Values outside the given samples will use the outmost sample
   * as a constant extrapolation.
   */
  InterpolateData2DSpline(const std::vector<double>& x,
                          const std::vector<double>& y,
                          const std::vector<double>& z);

  /// Destructor
  ~InterpolateData2DSpline();

  /**
   * Calculate bicubic interpolation for given x and y.
   *
   *  \param xi Interpolation argument in first dimension.
   *  \param yi Interpolation argument in second dimension.
   *  \return Interpolated value.
   */
  double operator()(double xi, double yi) const;

 private:
  /// First x value.
  double first_x_;
  /// Last x value.
  double last_x_;
  /// First y value.
  double first_y_;
  /// Last y value.
  double last_y_;

  /// GSL iterator for interpolation lookups in x direction.
  gsl_interp_accel* xacc_;
  /// GSL iterator for interpolation lookupin y direction.
  gsl_interp_accel* yacc_;
  /// GSL spline in 2D.
  gsl_spline2d* spline_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_INTERPOLATION2D_H_
