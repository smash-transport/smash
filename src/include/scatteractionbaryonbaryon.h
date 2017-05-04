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
  /** Determine the parametrized string excitation cross section
   * for a baryon-baryon collision. */
  float string_cross_section() const override;
  /* There is no resonance formation out of two baryons: Return empty list. */
  CollisionBranchList resonance_cross_sections() override {
    return CollisionBranchList();
  }
  /** Find all inelastic 2->2 processes for this reaction. */
  CollisionBranchList two_to_two_cross_sections() override;

 private:
  /**
  * Calculate cross sections for resonance absorption
  * (i.e. NR->NN and ΔR->NN).
  *
  * \param[in] is_anti_particles Whether the colliding particles are antiparticles
  *
  * \return List of possible resonance absorption processes. Each element of the
  * list contains the types of the final-state particles and the cross section
  * for that particular process.
  */
  CollisionBranchList bar_bar_to_nuc_nuc(const bool is_anti_particles);

 protected:
  /**
   * Scattering matrix amplitude squared (divided by 16π) for resonance
   * production processes like NN → NR and NN → ΔR, where R is a baryon
   * resonance (Δ, N*, Δ*). Includes no spin or isospin factors.
   *
   * \param[in] srts sqrt(Mandelstam-s), i.e. collision CMS energy.
   * \param[in] type_a Type information for the first final-state particle.
   * \param[in] type_b Type information for the second final-state particle.
   * \param[in] twoI Twice the total isospin of the involved state.
   *
   * \return Matrix amplitude squared \f$ |\mathcal{M}(\sqrt{s})|^2/16\pi \f$.
   */
  static float nn_to_resonance_matrix_element(const double srts,
                                              const ParticleType &type_a,
                                              const ParticleType &type_b,
                                              const int twoI);

  /**
   * \ingroup logging
   * Writes information about this scatter action to the \p out stream.
   */
  void format_debug_output(std::ostream &out) const override;
};


}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONBARYONBARYON_H_
