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
  /**
   * Determine the elastic cross section for a nucleon-nucleon collision.
   * It is given by a parametrization of experimental data.
   *
   * \param[in] elast_par Elastic cross section parameter from the input file (not used here).
   *
   * \return A ProcessBranch object containing the cross section and
   * final-state IDs.
   */
  CollisionBranchPtr elastic_cross_section(float elast_par) override;
  /** Find all inelastic 2->2 processes for this reaction. */
  CollisionBranchList two_to_two_cross_sections() override;

 protected:
  /** Perform an elastic nucleon-nucleon scattering with
   * anisotropic angular distributions. */
  virtual void elastic_scattering() override;

  /**
   * Sample final state momenta (and masses) in an inelastic 2->2 collision,
   * possibly using anisotropic angular distributions.
   *
   * \throws InvalidResonanceFormation
   */
  virtual void sample_cms_momenta() override;

 private:
  /**
   * Calculate cross sections for single-resonance production from
   * nucleon-nucleon collisions (i.e. NN->NR).
   *
   * Checks are processed in the following order:
   * 1. Charge conservation
   * 2. Isospin factors (Clebsch-Gordan)
   * 3. Enough energy for all decay channels to be available for the resonance
   *
   * \param[in] type_particle1 Type information of the first incoming nucleon.
   * \param[in] type_particle2 Type information of the second incoming nucleon.
   *
   * \return List of resonance production processes possible in the collision
   * of the two nucleons. Each element in the list contains the type(s) of the
   * final state particle(s) and the cross section for that particular process.
   */
  CollisionBranchList nuc_nuc_to_nuc_res(const ParticleType &type_particle1,
                                       const ParticleType &type_particle2);
};


}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONNUCLEONNUCLEON_H_
