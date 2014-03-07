/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_RANDOM_H_
#define SRC_INCLUDE_RANDOM_H_

#include<random>

namespace Smash {


/**
 * \todo{Random number generator is a global object at the moment. This
 * needs to be revised for multi-threading}
 *
 **/
class rng_dist {
 public:
  rng_dist() 
     : engine_(time(nullptr)),
       cos_like_(-1.0, 1.0),
       phi_like_(0.0, 2 * M_PI),
       uniform_(0.0, 1.0),
       // canonical_(0.0, 1.0),
       exponential_(1.0) {
  }
  rng_dist(int use_seed) 
     : engine_(use_seed),
       cos_like_(-1.0, 1.0),
       phi_like_(0.0, 2 * M_PI),
       uniform_(0.0, 1.0),
       // canonical_(0.0, 1.0),
       exponential_(1.0) {
  }

  void seed(int use_seed) {
    engine_.seed(use_seed);
  }
  double phi_like() {
    return phi_like_(engine_);
  }
  double cos_like() {
    return cos_like_(engine_);
  }
  double uniform() {
    return uniform_(engine_);
  }
  double uniform(const double& a, const double& b) {
    set_uniform(a,b);
    return uniform();
  }
  double canonical() {
    return std::generate_canonical<double,52>(engine_);
  }
  void set_uniform(const double& a, const double& b) {
    uniform_.param(std::uniform_real_distribution<>::param_type(a,b));
  }
  double exponential() {
    return exponential_(engine_);
  }

 private:
  /** The random number engine used in SMASH. **/
  std::ranlux48 engine_;
  /** provides cosine-like random numbers (uniform distribution in
   * [-1..1) )
   **/
  std::uniform_real_distribution<double> cos_like_;
  /** provides phi-like random numbers (uniform distribution in
   * [0..2*pi) )
   **/
  std::uniform_real_distribution<double> phi_like_;
  /** provides uniform random numbers **/
  std::uniform_real_distribution<double> uniform_;
  /** provides uniform random numbers between 0.0 and 1.0
   **/
  //std::generate_canonical<double,52> canonical_;
  //std::uniform_real_distribution<double> canonical_;
  /** provides exponential random numbers (\f$P(t) = \exp{-t}\f$) **/
  std::uniform_real_distribution<double> exponential_;
};

extern rng_dist rng;

}  // namespace Smash

#endif  // SRC_INCLUDE_RANDOM_H_
