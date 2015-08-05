/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_TABULATION_H_
#define SRC_INCLUDE_TABULATION_H_

#include <vector>
#include <map>
#include <memory>

#include "forwarddeclarations.h"
#include "particletype.h"

namespace Smash {

/**
 * A class for storing a one-dimensional lookup table of floating-point values.
 */
class Tabulation {
 public:
  /** Construct a new tabulation object.
   * \param x_min lower bound of tabulation domain
   * \param range range (x_max-x_min) of tabulation domain
   * \param num number of intervals (the number of tabulated points is actually num+1)
   * \param f one-dimensional function f(x) which is supposed to be tabulated
   */
  Tabulation(float x_min, float range, int num,
             std::function<double(float)> f);
  /** Look up a value from the tabulation (without any interpolation, simply
   * using the closest tabulated value). If x is below the lower tabulation
   * bound we return 0, if it is above the upper bound we return the tabulated
   * value at the upper bound. */
  float get_value_step(float x) const;
  /** Look up a value from the tabulation using linear interpolation.
   * If x is above the upper bound, we use linear extrapolation of the two
   * highest tabulated points. */
  float get_value_linear(float x) const;

 protected:
  // vector for storing tabulated values
  std::vector<float> values_;
  // lower bound and inverse step size 1/dx for tabulation
  const float x_min_, inv_dx_;
};


// a map for storing the tabulation used for the NN->NR cross sections
static std::map<int, TabulationPtr> XS_tabulation;


}  // namespace Smash

#endif  // SRC_INCLUDE_TABULATION_H_
