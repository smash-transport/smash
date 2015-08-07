/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ADAPTIVEREJECTIONSAMPLER_H_
#define SRC_INCLUDE_ADAPTIVEREJECTIONSAMPLER_H_

#include <vector>
#include <list>
#include <functional>
#include "random.h"

namespace Smash {

namespace Rejection {
/*x, expy=f(x), y=log(f(x)) coordinates used in AdaptiveRejectionSampler*/
typedef struct point {
    float x, y;
    float expy;
}Point;

std::ostream &operator<<(std::ostream &out, const Point &p);

/*lines used to define the upper bounds and lower bounds*/
typedef struct line {
    float m, b;                // f(x) = m*x+b
    float eval(float x) {
        return m*x + b;
    }
}Line;

/** Envelope to hold one piece of upper bounds*/
typedef struct envelope {
    Point left_point, right_point;
    Line  piecewise_linear_line;
}Envelope;

std::ostream &operator<<(std::ostream &out, const Line &l);

/**
 *Adaptive Rjection Sampling used for thermal, juttner,
 *bose-einstein, fermi-dirac and woods-saxon
 *distributions. They are all log concave distributions.
 *<a href="https://en.wikipedia.org/wiki/Rejection_sampling">see wikipedia</a>.
 *Here is an example of AdaptiveRejectionSampler usage:
 *\code
 *
 *float woods_saxon_dist(float r, float radius, float diffusion)
 *{
 *    return r*r/(exp((r-radius)/diffusion)+1.0);
 *}
 *
 *int main() {
 *    using namespace rejection;
 *    float radius = 6.4;
 *    float diffusion = 0.54;
 *    float xmin = 0.0;
 *    float xmax = 15.0;
 *    AdaptiveRejectionSampler sampler([&](float x) {
 *        return woods_saxon_dist(x, radius, diffusion);}
 *        ,xmin, xmax);
 *
 *    float x = sampler.get_one_sample();
 *    return 0;
 *}
 *\endcode
 */

class AdaptiveRejectionSampler {
 public:
     /* distribution function f(x) for sampling
      * (arguments are hiden by lambda functions)
      */
     std::function<float(float)> f_;

     /* The left end of the range */
     float xmin_ = 0.0;

     /* The right end of the range */
     float xmax_ = 15.0;

     /* Maximum refine loops to avoid further adaptive updates.*/
     int max_refine_loops_ = 40;

     /* Num of points to initialize the upper bound*/
     int init_npoint_ = 10;

     /* points on distribution function curve with (x,logf(x),f(x))*/
     std::vector<Point> points_;

     /* intersections between each pair of neighboring upper bounds*/
     std::vector<Point> inters_;

     /* scants that connect all the points on distribution function*/
     std::vector<Line>  scants_;

     /* store the upper bounds to make get_one_sample faster */
     std::vector<Envelope> upper_bounds_;

     /* areas below each piece of upper bound*/
     std::vector<float> areas_;

     /* areas_ list as weight for the discrete_distribution_ */
     Random::discrete_dist<float> discrete_distribution_;

     AdaptiveRejectionSampler() = default;

     /**
      *Constructor for adaptive rejection sampling
      *\param func: function pointer for the distribution function
      *\param xmin: minimum x in sampling f(x)
      *\param xmax: maximum x in sampling f(x)
      */
     AdaptiveRejectionSampler(std::function<float(float)> func,
             float xmin, float xmax);

     /*reset max refine loops for AdaptiveRejectionSampler*/
     void reset_max_refine_loops(const int new_max_refine_loops);

     /*sample one x from distribution function f(x) */
     float get_one_sample();

 private:
     /*initialize scants with 10 points between xmin and xmax by
      * default or with user provided xlist
      * */
     void init_scant();

     /*initialize intersections with initial scants*/
     void init_inter();

     // get area below piecewise exponential upper bound
     void update_area();

     // sample j from discrete distribution with weight by area list
     int sample_j();

     // return the upper bound of the j'th piece of area
     inline Line upper(int j);

     // return the lower bound of the j'th piece of area
     inline Line lower(int j);

     // sample x in the j'th piece
     float sample_x(int j);

     // r<exp(lower-upper)
     inline bool squeezing_test(const float & x, const int & j,
             const float & rand);

     // r<func/exp(upper)
     inline bool rejection_test(const float & x, const int & j,
             const float & rand);

     /// calc line from 2 points
     inline Line create_line(Point p0, Point p1);

     /// calc intersection from 2 scant lines
     inline Point create_inter(Line l0, Line l2);

     /// calc left and right most points in upper bounds
     inline void create_leftend();
     inline void create_rightend();

     // refine scants, intersections, after one rejection
     void adaptive_update(const int j, const Point & new_rejection);
};


}  // end namespace Rejection

}  // end namespace Smash

#endif  // SRC_INCLUDE_ADAPTIVEREJECTIONSAMPLER_H_
