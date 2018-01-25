
/*
 *
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionbaryonbaryon.h"

#include "include/clebschgordan.h"
#include "include/constants.h"
#include "include/cxx14compat.h"
#include "include/parametrizations.h"

namespace smash {

double ScatterActionBaryonBaryon::high_energy_cross_section() const {
  const PdgCode &pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode &pdg_b = incoming_particles_[1].type().pdgcode();
  const double s = mandelstam_s();

  /* Currently all BB collisions use the nucleon-nucleon parametrizations. */
  if (pdg_a == pdg_b) {
    return pp_high_energy(s);  // pp, nn
  } else if (pdg_a.is_antiparticle_of(pdg_b)) {
    return ppbar_high_energy(s);  // ppbar, nnbar
  } else if (pdg_a.antiparticle_sign() * pdg_b.antiparticle_sign() == 1) {
    return np_high_energy(s);  // np, nbarpbar
  } else {
    return npbar_high_energy(s);  // npbar, nbarp
  }
}

double ScatterActionBaryonBaryon::string_hard_cross_section() const {
  // const PdgCode &pdg_a = incoming_particles_[0].type().pdgcode();
  // const PdgCode &pdg_b = incoming_particles_[1].type().pdgcode();
  const double s = mandelstam_s();

  /**
   * Currently nucleon-nucleon cross section is used for all case.
   * This will be changed later by applying additive quark model.
   */
  return NN_string_hard(s);
}

void ScatterActionBaryonBaryon::format_debug_output(std::ostream &out) const {
  out << "Baryon-Baryon ";
  ScatterAction::format_debug_output(out);
}

}  // namespace smash
