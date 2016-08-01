/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONNUCLEONKAON_H_
#define SRC_INCLUDE_SCATTERACTIONNUCLEONKAON_H_

#include "scatteraction.h"
#include "scatteractionbaryonmeson.h"

namespace Smash {


/**
 * \ingroup action
 * ScatterActionNucleonKaon is a special ScatterActionBaryonMeson which represents the
 * scattering of a nucleon and a kaon.
 */
class ScatterActionNucleonKaon : public ScatterActionBaryonMeson {
 public:
  /* Inherit constructor. */
  using ScatterActionBaryonMeson::ScatterActionBaryonMeson;
  /**
   * Determine the elastic cross section for a nucleon-kaon collision.
   * It is given by a parametrization of experimental data.
   */
  float elastic_parametrization() override;
  /** Find all inelastic 2->2 processes for this reaction. */
  CollisionBranchList two_to_two_cross_sections() override;

 protected:
  /**
   * \ingroup logging
   * Writes information about this scatter action to the \p out stream.
   */
  void format_debug_output(std::ostream &out) const override;

 private:
  /**
   * Calculate cross sections for strangeness exchange in nucleon-kaon
   * collisions.
   */
  CollisionBranchList two_to_two_inel(const ParticleType &type_particle_a,
                                      const ParticleType &type_particle_b);
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONNUCLEONKAON_H_
