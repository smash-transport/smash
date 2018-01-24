/*
 *
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionnucleonpion.h"

#include "include/clebschgordan.h"
#include "include/cxx14compat.h"
#include "include/parametrizations.h"
#include "include/pdgcode_constants.h"

namespace smash {

void ScatterActionNucleonPion::format_debug_output(std::ostream &out) const {
  out << "Nucleon-Pion  ";
  ScatterAction::format_debug_output(out);
}

}  // namespace smash
