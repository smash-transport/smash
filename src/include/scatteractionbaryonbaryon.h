/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONBARYONBARYON_H_
#define SRC_INCLUDE_SCATTERACTIONBARYONBARYON_H_

#include "scatteraction.h"

namespace Smash {


/**
 * \ingroup action
 * ScatterActionBaryonBaryon is a special ScatterAction which represents the
 * scattering of two baryons.
 */
class ScatterActionBaryonBaryon : public ScatterAction {
 public:
  /* Inherit constructor. */
  using ScatterAction::ScatterAction;
  /** Determine the parametrized total cross section
   * for a baryon-baryon collision. */
  virtual float total_cross_section() const override;
  /* There is no resonance formation out of two baryons: Return empty list. */
  CollisionBranchList resonance_cross_sections() override {
    return CollisionBranchList();
  }
  /** Find all inelastic 2->2 processes for this reaction. */
  CollisionBranchList two_to_two_cross_sections() override;

 private:
  /**
  * Calculate cross sections for resonance absorption on a nucleon
  * (i.e. NR->NN).
  *
  * \param[in] type_particle1 Type information of the first incoming nucleon.
  * \param[in] type_particle2 Type information of the second incoming nucleon.
  *
  * \return List of resonance absorption processes possible in the collision
  * with a nucleon. Each element in the list contains the type(s) of the
  * final state particle(s) and the cross section for that particular process.
  */
  CollisionBranchList nuc_res_to_nuc_nuc(const ParticleType &type_particle1,
                                       const ParticleType &type_particle2);

 protected:
  /**
   * \ingroup logging
   * Writes information about this scatter action to the \p out stream.
   */
  void format_debug_output(std::ostream &out) const override;
};


}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONBARYONBARYON_H_
