/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_TABULATION_H_
#define SRC_INCLUDE_TABULATION_H_

#include <vector>

#include "particletype.h"

namespace Smash {

/** Parameters for GSL integration. */
struct IntegParam {
  const ParticleType &type;  // type of daughter resonance
  double m2;                 // mass of stable particle
  double srts;               // sqrt(s) = mass of decaying resonance
  int L;                     // angular momentum
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
   * \param N number of tabulation points
   * \param ip integration parameters
   * \param f integrand function
   */
  Tabulation(float x_min, float range, unsigned int N,
             IntegParam ip, IntegrandFunction f);
  /// Look up a value from the tabulation.
  float get_value(float x) const;

 protected:
  /** Calculate a value to be stored in the table
   * (by numerically solving an integral). */
  float calculate_value(float x);
  std::vector<float> values_;  // vector for storing tabulated values
  float x_min_, dx_;           // minimum mass and step size for tabulation
  IntegParam ip_;              // integration parameters
  IntegrandFunction func_;     // integrand function
};

}  // namespace Smash

#endif  // SRC_INCLUDE_TABULATION_H_
