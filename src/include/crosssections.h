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

/** Contains the parametrizations of Cross-Sections
 *
 * This Class is deprecated and will be replaced in DecayAction
 */
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
  /** Returns the elastic cross-section of a collision of the particles
   * \a a and \a b.
   *
   * \param data_a Data of first particle (\a a)
   * \param data_b Data of second particle (\a b)
   */
  float elastic(const ParticleData &data_a, const ParticleData &data_b);

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
  /** Computes the kinematic variables used in the further calculations.
   *
   * \param data_a Data of first particle (\a a)
   * \param data_b Data ID of second particle (\a b)
   *
   * This function sets
   * \see squared_mass_a_
   * \see squared_mass_b_
   * \see mandelstam_s_
   * \see p_lab_
   */
  void compute_kinematics(const ParticleData &data_a,
                          const ParticleData &data_b);
  /** Returns the total (elastic + inelastic) cross-section of a
   * collision of the particles \a a and \a b.
   *
   * \param pdg_a PDG code of first particle (\a a)
   * \param pdg_b PDG code of second particle (\a b)
   */
  float total(const PdgCode &pdg_a, const PdgCode &pdg_b) const;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_CROSSSECTIONS_H_
