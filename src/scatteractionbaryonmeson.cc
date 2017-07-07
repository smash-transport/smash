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

float ScatterActionBaryonMeson::high_energy_cross_section() const {
  const PdgCode &pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode &pdg_b = incoming_particles_[1].type().pdgcode();
  const double s = mandelstam_s();

  /* Currently only include pion nucleon interaction. */
  if ((pdg_a == pdg::pi_p && pdg_b == pdg::p)
      ||(pdg_b == pdg::pi_p && pdg_a == pdg::p)
      || (pdg_a == pdg::pi_m && pdg_b == pdg::n)
      || (pdg_b == pdg::pi_m && pdg_a == pdg::n)) {
    return piplusp_high_energy(s);     // pi+ p, pi- n
  } else if ((pdg_a == pdg::pi_m && pdg_b == pdg::p)
            || (pdg_b == pdg::pi_m && pdg_a == pdg::p)
            || (pdg_a == pdg::pi_p && pdg_b == pdg::n)
            || (pdg_b == pdg::pi_p && pdg_a == pdg::n)) {
    return piminusp_high_energy(s);   // pi- p, pi+ n
  } else {
     return 0;
  }
}

}  // namespace Smash
