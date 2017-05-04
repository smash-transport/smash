/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionbaryonmeson.h"
#include "include/parametrizations.h"

namespace Smash {


void ScatterActionBaryonMeson::format_debug_output(std::ostream &out) const {
  out << "Baryon-Meson  ";
  ScatterAction::format_debug_output(out);
}

float ScatterActionBaryonMeson::string_cross_section() const {
  const PdgCode &pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode &pdg_b = incoming_particles_[1].type().pdgcode();
  const double s = mandelstam_s();

  /* Currently only include pion nucleon interaction. */
  if ((pdg_a == 211 && pdg_b == 2212) ||( pdg_b == 211 && pdg_a == 2212)
      || (pdg_a == -211 && pdg_b == 2112) || (pdg_b == -211 && pdg_a == 2112)) {
    return piplusp_string(s);     // pi+ p, pi- n
  } else if ((pdg_a == -211 && pdg_b == 2212) ||( pdg_b == -211 && pdg_a == 2212)
      || (pdg_a == 211 && pdg_b == 2112) || (pdg_b == 211 && pdg_a == 2112)) {
    return piminusp_string(s);   // pi- p, pi+ n
  } else {
     return 0;
  }
}

}  // namespace Smash
