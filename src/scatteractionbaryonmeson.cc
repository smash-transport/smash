/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionbaryonmeson.h"

namespace Smash {


void ScatterActionBaryonMeson::format_debug_output(std::ostream &out) const {
  out << "Baryon-Meson  ";
  ScatterAction::format_debug_output(out);
}


}  // namespace Smash
