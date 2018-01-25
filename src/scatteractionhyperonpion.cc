/*
 *
 *    Copyright (c) 2016-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionhyperonpion.h"

#include "include/cxx14compat.h"
#include "include/parametrizations.h"
#include "include/pdgcode_constants.h"

namespace smash {

void ScatterActionHyperonPion::format_debug_output(std::ostream& out) const {
  out << "Hyperon-Pion  ";
  ScatterAction::format_debug_output(out);
}


}  // namespace smash
