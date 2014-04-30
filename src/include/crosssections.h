/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_CROSSSECTIONS_H_
#define SRC_INCLUDE_CROSSSECTIONS_H_

#include "include/particles.h"

namespace Smash {

class CrossSections {
 public:
  /* Default constructor */
  explicit CrossSections(float elastic_parameter)
      : elastic_parameter_(elastic_parameter) {
  }
  /* Compute kinematics */
  void compute_kinematics(Particles *particles, int id_a, int id_b);
  /* Return cross sections */
  float elastic(Particles *particles, int id_a, int id_b) const;
  float total(Particles *particles, int id_a, int id_b) const;

  /// Resets the parameters to the values that were set in the constructor.
  void reset() {
    squared_mass_a_ = -1.;
    squared_mass_b_ = -1.;
    mandelstam_s_ = -1.;
    p_lab_ = -1.;
  }

 private:
  /* Elastic cross section parameter */
  const float elastic_parameter_ = 0.0;
  /* Mass of the first particle */
  float squared_mass_a_ = -1.0;
  /* Mass of the second particle */
  float squared_mass_b_ = -1.0;
  /* Mandelstam s of the collision (= CMS energy squared) */
  double mandelstam_s_ = -1.0;
  /* "Beam" momentum */
  double p_lab_ = -1.0;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_CROSSSECTIONS_H_
