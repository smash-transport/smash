
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
