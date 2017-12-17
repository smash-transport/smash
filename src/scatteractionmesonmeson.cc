/*
 *
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionmesonmeson.h"
#include "include/parametrizations.h"

namespace smash {

void ScatterActionMesonMeson::format_debug_output(std::ostream &out) const {
  out << " Meson-Meson  ";
  ScatterAction::format_debug_output(out);
}

double ScatterActionMesonMeson::string_hard_cross_section() const {
  // const PdgCode &pdg_a = incoming_particles_[0].type().pdgcode();
  // const PdgCode &pdg_b = incoming_particles_[1].type().pdgcode();
  const double s = mandelstam_s();

  /**
   * Currently pion-pion cross section is used for all case.
   * This will be changed later by applying additive quark model.
   */
  return pipi_string_hard(s);
}

}  // namespace smash
