/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONNUCLEONNUCLEON_H_
#define SRC_INCLUDE_SCATTERACTIONNUCLEONNUCLEON_H_

#include <utility>

#include "scatteractionbaryonbaryon.h"

namespace Smash {


/**
 * \ingroup action
 * ScatterActionNucleonNucleon is a special ScatterAction which represents the
 * scattering of two nucleons.
 */
class ScatterActionNucleonNucleon : public ScatterActionBaryonBaryon {
 public:
  /* Inherit constructor. */
  using ScatterActionBaryonBaryon::ScatterActionBaryonBaryon;
  /** Determine the (parametrized) elastic cross section for a
   * nucleon-nucleon collision. */
  double elastic_parametrization() override;
  /** Find all inelastic 2->2 processes for this reaction.
   * Calculate cross sections for resonance production from
   * nucleon-nucleon collisions (i.e. N N -> N R, N N -> Delta R).
   *
   * Checks are processed in the following order:
   * 1. Charge conservation
   * 2. Isospin factors (Clebsch-Gordan)
   * 3. Enough energy for all decay channels to be available for the resonance
   *
   * \return List of resonance production processes possible in the collision
   * of the two nucleons. Each element in the list contains the type(s) of the
   * final state particle(s) and the cross section for that particular process.
   */
  CollisionBranchList two_to_two_cross_sections() override;

 protected:
  /**
   * Sample final-state angles in a 2->2 collision (possibly anisotropic).
   *
   * \throws InvalidResonanceFormation
   */
  void sample_angles(std::pair<double, double> masses) override;

 private:
  /**
   * Utility function to avoid code replication in two_to_two_cross_sections
   */
  template<class IntegrationMethod>
  CollisionBranchList find_xsection_from_type(
                                         const ParticleTypePtrList &type_res_1,
                                         const ParticleTypePtrList &type_res_2,
                                         const IntegrationMethod integrator);
};


}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONNUCLEONNUCLEON_H_
