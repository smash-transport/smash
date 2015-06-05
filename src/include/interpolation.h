/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_INTERPOLATION_H_
#define SRC_INCLUDE_INTERPOLATION_H_

#include <vector>

class InterpolateLinear {
  public:
    double m, n;
    /// Linear interpolation given two points (x0, y0) and (x1, y1).
    ///
    /// Returns the interpolation function.
    InterpolateLinear(double x0, double y0, double x1, double y1);
    /// Calculate linear interpolation at x.
    double operator()(double x) const;
};

class InterpolateData {
  public:
    /// Interpolate function f given discrete samples f(x_i) = y_i.
    //
    /// Returns the interpolation function.
    ///
    /// Piecewise linear interpolation is used.
    /// Values outside the given samples will use the outmost linear interpolation.
    InterpolateData(const std::vector<double>& x, const std::vector<double>& y);
    /// Calculate interpolation of f at x.
    double operator()(double x) const;
  private:
    std::vector<double> x;
    std::vector<InterpolateLinear> f;
};

#endif  // SRC_INCLUDE_INTERPOLATION_H_
