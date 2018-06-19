/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ADAPTIVEREJECTIONSAMPLER_H_
#define SRC_INCLUDE_ADAPTIVEREJECTIONSAMPLER_H_

#include <functional>
#include <list>
#include <vector>

#include "random.h"

namespace smash {

namespace rejection {
/**
 * x, expy=f(x), y=log(f(x)) coordinates used in AdaptiveRejectionSampler
 */
struct Point {
  double x;
  double y;
  double expy;
};

std::ostream &operator<<(std::ostream &out, const Point &p);

/**
 * Lines used to define the upper bounds and lower bounds
 *
 * f(x) = m*x + b
 */
struct Line {
  double m;
  double b;
  double eval(double x) { return m * x + b; }
};

/**
 * Envelope to hold one piece of upper bounds
 */
struct Envelope {
  Point left_point;
  Point right_point;
  Line piecewise_linear_line;
};

std::ostream &operator<<(std::ostream &out, const Line &l);

/**
 * Adaptive Rjection Sampling used for thermal, Juttner, Bose-Einstein,
 * Fermi-Dirac and Woods-Saxon distributions. They are all log concave
 * distributions.
 * <a href="https://en.wikipedia.org/wiki/Rejection_sampling">See wikipedia</a>.
 *
 * Here is an example of AdaptiveRejectionSampler usage:
 * \code
 *
 * double woods_saxon_dist(double r, double radius, double diffusion)
 * {
 *     return r*r/(exp((r-radius)/diffusion)+1.0);
 * }
 *
 * int main() {
 *     using namespace rejection;
 *     double radius = 6.4;
 *     double diffusion = 0.54;
 *     double xmin = 0.0;
 *     double xmax = 15.0;
 *     AdaptiveRejectionSampler sampler(
 *         [&](double x) {
 *             return woods_saxon_dist(x, radius, diffusion);
 *         }, xmin, xmax);
 *
 *     double x = sampler.get_one_sample();
 *     return 0;
 * }
 * \endcode
 */
/// \todo unused
class AdaptiveRejectionSampler {
 public:
  /// Distribution function f(x) for sampling.
  std::function<double(double)> f_;

  /// The left end of the range
  double xmin_ = 0.0;

  /// The right end of the range.
  double xmax_ = 15.0;

  /// Maximum refine loops to avoid further adaptive updates.
  int max_refine_loops_ = 40;

  /// Num of points to initialize the upper bound.
  int init_npoint_ = 10;

  /// Points on distribution function curve with (x, logf(x), f(x)).
  std::vector<Point> points_;

  /// Intersections between each pair of neighboring upper bounds.
  std::vector<Point> inters_;

  /// Scants that connect all the points on distribution function.
  std::vector<Line> scants_;

  /// Store the upper bounds to make get_one_sample faster.
  std::vector<Envelope> upper_bounds_;

  /// Areas below each piece of upper bound.
  std::vector<double> areas_;

  /// areas_ list as weight for the discrete_distribution_.
  random::discrete_dist<double> discrete_distribution_;

  AdaptiveRejectionSampler() = default;

  /**
   * Constructor for adaptive rejection sampling.
   *
   * \param func: function pointer for the distribution function
   * \param xmin: minimum x in sampling f(x)
   * \param xmax: maximum x in sampling f(x)
   */
  AdaptiveRejectionSampler(std::function<double(double)> func, double xmin,
                           double xmax);

  // Reset max refine loops.
  // void reset_max_refine_loops(const int new_max_refine_loops);

  /// Sample one x from distribution function f(x).
  double get_one_sample();

 private:
  /**
   * Initialize scants with 10 points between xmin and xmax by
   * default or with user provided xlist.
   */
  void init_scant();

  /// Initialize intersections with initial scants.
  void init_inter();

  /// Get area below piecewise exponential upper bound.
  void update_area();

  /// Sample j from discrete distribution with weight by area list.
  int sample_j();

  /// Return the upper bound of the j'th piece of area.
  inline Line upper(int j);

  /// Return the lower bound of the j'th piece of area.
  inline Line lower(int j);

  /// Sample x in the j'th piece.
  double sample_x(int j);

  /// r < exp(lower-upper)
  inline bool squeezing_test(const double &x, const int &j, const double &rand);

  /// r < func/exp(upper)
  inline bool rejection_test(const double &x, const int &j, const double &rand);

  /// Calculate line from 2 points.
  inline Line create_line(Point p0, Point p1);

  /// Calculate intersection from 2 scant lines.
  inline Point create_inter(Line l0, Line l2);

  /// Calculate left-most points in upper bounds.
  inline void create_leftend();

  /// Calculate right-most points in upper bounds.
  inline void create_rightend();

  /// Refine scants, intersections, after one rejection.
  void adaptive_update(const int j, const Point &new_rejection);
};

}  // end namespace rejection

}  // end namespace smash

#endif  // SRC_INCLUDE_ADAPTIVEREJECTIONSAMPLER_H_
