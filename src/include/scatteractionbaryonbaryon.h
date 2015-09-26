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
  * Calculate cross sections for resonance absorption
  * (i.e. NR->NN and DeltaR->NN).
  *
  * \param[in] type_a Type information of the first incoming baryon.
  * \param[in] type_b Type information of the second incoming baryon.
  *
  * \return List of possible resonance absorption processes. Each element of the
  * list contains the types of the final-state particles and the cross section
  * for that particular process.
  */
  CollisionBranchList bar_bar_to_nuc_nuc(const ParticleType &type_a,
                                         const ParticleType &type_b);

 protected:
  /**
   * Scattering matrix amplitude squared for resonance production processes like
   * \f$ NN \rightarrow NR \f$  and \f$ NN \rightarrow \Delta R \f$,
   * where R is a baryon resonance (\f$ \Delta, N^*, \Delta^* \f$).
   * Includes a spin factor \f$ (2S_a+1)(2S_b+1) \f$, but no isospin factors.
   *
   * \param[in] srts sqrt(Mandelstam-s), i.e. collision CMS energy.
   * \param[in] type_a Type information for the first final-state particle.
   * \param[in] type_b Type information for the second final-state particle.
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
