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

/** Parameters for GSL integration. */
struct IntegParam {
  const ParticleTypePtr type;  // type of final-state resonance
  double m2;                   // mass of stable final-state particle
  double srts;                 // sqrt(s)
  int L;                       // angular momentum
};

typedef double (*IntegrandFunction)(double, void*);

/**
 * A class for storing a one-dimensional lookup table of floating-point values,
 * which are obtained by integration.
 */
class Tabulation {
 public:
  /** Construct a new tabulation object.
   * \param x_min lower bound of tabulation domain
   * \param range range (x_max-x_min) of tabulation domain
   * \param num_points number of tabulation points
   * \param ip integration parameters
   * \param f integrand function
   */
  Tabulation(float x_min, float range, int num_points,
             IntegParam ip, IntegrandFunction f);
  /** Look up a value from the tabulation. If x is below the lower tabulation
   * bound we return 0, if it is above the upper bound we return the tabulated
   * value at the upper bound. */
  float get_value(float x) const;

 protected:
  /** Calculate a value to be stored in the table
   * (by numerically solving an integral). */
  float calculate_value(float x, IntegParam ip, IntegrandFunction f);
  std::vector<float> values_;   // vector for storing tabulated values
  const float x_min_, inv_dx_;  // lower bound and inverse step size 1/dx for tabulation
};


// a map for storing the integral tabulation used for the NN->NR cross sections
static std::map<int, TabulationPtr> XS_tabulation;


}  // namespace Smash

#endif  // SRC_INCLUDE_TABULATION_H_
