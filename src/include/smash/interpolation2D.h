/*
 *
 *    Copyright (c) 2020,2022,2026
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_INTERPOLATION2D_H_
#define SRC_INCLUDE_SMASH_INTERPOLATION2D_H_

#include <vector>

#include "gsl/gsl_spline2d.h"

#include "smash/forwarddeclarations.h"

namespace smash {

/// Represent a bicubic spline interpolation.
class InterpolateData2DSpline {
 public:
  /**
   * Interpolate function f given discrete samples f(x_i, y_i) = z_i.
   * A bicubic spline interpolation is used.
   *
   * \param x x-values.
   * \param y y-values.
   * \param z z-values
   * \param extrapolation_type Type of extrapolation for requested x_i and y_i
   *                           values that are out of bounds. Extrapolation is
   *                           by default disabled. Possible types are
   *                           <tt>None</tt> and <tt>Constant</tt>.
   *
   * \return The interpolation function.
   * \throw std::out_of_range if values outside of the boundaries of the
   *                          underlying data are tried to be accessed and
   *                          extrapolation is disabled.
   * \throw std::invalid_argument if unsupported extrapolation type is
   *                              requested.
   *
   */
  InterpolateData2DSpline(
      const std::vector<double>& x, const std::vector<double>& y,
      const std::vector<double>& z,
      const ExtrapolationType extrapolation_type = ExtrapolationType::None);

  /// Destructor
  ~InterpolateData2DSpline();

  /**
   * Calculate bicubic interpolation for given x and y.
   *
   * \param xi Interpolation argument in first dimension.
   * \param yi Interpolation argument in second dimension.
   *
   * \return Interpolated value.
   * \throw std::out_of_range if values outside of the boundaries of the
   *                          underlying data are tried to be accessed and
   *                          extrapolation is disabled.
   */
  double operator()(double xi, double yi) const;

 private:
  /// Extrapolation type.
  ExtrapolationType extrapolation_type_;
  /// First x value of underlying data.
  double first_x_;
  /// Last x value of underlying data.
  double last_x_;
  /// First y value of underlying data.
  double first_y_;
  /// Last y value of underlying data.
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
