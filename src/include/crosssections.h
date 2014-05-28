/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_CROSSSECTIONS_H_
#define SRC_INCLUDE_CROSSSECTIONS_H_

#include "forwarddeclarations.h"

namespace Smash {

class CrossSections {
 public:
  /** Constructs a new cross-sections object
   *
   * \param elastic_parameter sets the fall-back elastic cross-section,
   * \see elastic_parameter_
   */
  explicit CrossSections(float elastic_parameter)
      : elastic_parameter_(elastic_parameter) {
  }
  /** computes the kinematic variables used in the further calculations
   *
   * \param particles Particle list
   * \param id_a Unique ID of first particle (\a a)
   * \param id_b Unique ID of second particle (\a b)
   *
   * This function sets
   * \see squared_mass_a_
   * \see squared_mass_b_
   * \see mandelstam_s_
   * \see p_lab_
   */
  void compute_kinematics(Particles *particles, int id_a, int id_b);
  /** returns the elastic cross-section of a collision of the particles
   * \a a and \a b.
   *
   * \param particles Particle list
   * \param id_a Unique ID of first particle (\a a)
   * \param id_b Unique ID of second particle (\a b)
   */
  float elastic(Particles *particles, int id_a, int id_b) const;
  /** returns the total (elastic + inelastic) cross-section of a
   * collision of the particles \a a and \a b.
   *
   * \param particles Particle list
   * \param id_a Unique ID of first particle (\a a)
   * \param id_b Unique ID of second particle (\a b)
   */
  float total(Particles *particles, int id_a, int id_b) const;

  /// Resets the parameters to the default values.
  void reset() {
    squared_mass_a_ = -1.;
    squared_mass_b_ = -1.;
    mandelstam_s_ = -1.;
    p_lab_ = -1.;
  }

 private:
  /** Elastic Cross section, in mb
   *
   * This is the cross-section that is returned for elastic collisions.
   *
   */
  const float elastic_parameter_ = 0.0;
  /// mass of first particle (\a a), squared
  float squared_mass_a_ = -1.0;
  /// mass of second particle (\a b), squared
  float squared_mass_b_ = -1.0;
  /** Mandelstam s of the collision (= total center-of-mass energy
   * squared 
   */
  double mandelstam_s_ = -1.0;
  /// Momentum of particle \a a in center-of-mass-frame
  double p_lab_ = -1.0;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_CROSSSECTIONS_H_
