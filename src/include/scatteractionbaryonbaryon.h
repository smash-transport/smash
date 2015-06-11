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
  float total_cross_section() const override;
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
   * Scattering matrix amplitude squared for \f$ NN \rightarrow NR \f$ processes,
   * where R is a baryon resonance (\f$ \Delta, N^*, \Delta^* \f$).
   * Includes a spin factor \f$ (2S_a+1)(2S_b+1) \f$, but no isospin factors.
   *
   * \param[in] srts sqrt(Mandelstam-s), i.e. collision CMS energy.
   * \param[in] type_a Type information for the first final state particle.
   * \param[in] type_b Type information for the second final state particle.
   *
   * \return Matrix amplitude squared \f$ |\mathcal{M}(\sqrt{s})|^2/16\pi \f$.
   */
  float nn_to_resonance_matrix_element(const double srts,
                                       const ParticleType &type_a,
                                       const ParticleType &type_b) const;

  /**
   * \ingroup logging
   * Writes information about this scatter action to the \p out stream.
   */
  void format_debug_output(std::ostream &out) const override;
};


}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONBARYONBARYON_H_
