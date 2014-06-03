/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/decaymodes.h"

#include "include/constants.h"
#include "include/pdgcode.h"
#include "include/processbranch.h"

#include <cstdio>

namespace Smash {

void DecayModes::add_mode(std::vector<PdgCode> pdg_list, float ratio) {
  if (pdg_list.size() < 2) {
    throw InvalidDecay(
        "DecayModes::add_mode was instructed to add a decay mode with less "
        "than 2 branches. This is an invalid input.");
  }
  ProcessBranch branch;
  branch.set_particles(std::move(pdg_list));
  branch.set_weight(ratio);
  decay_modes_.push_back(branch);
}

void DecayModes::renormalize(float renormalization_constant) {
  if (renormalization_constant < really_small) {
    printf("Warning: Extremely small renormalization constant: %g\n",
           renormalization_constant);
    printf("Skipping the renormalization.\n");
  } else {
    printf("Renormalizing decay modes with %g \n", renormalization_constant);
    float new_sum = 0.0;
    for (std::vector<ProcessBranch>::iterator mode
           = decay_modes_.begin(); mode != decay_modes_.end(); ++mode) {
      mode->set_weight(mode->weight() / renormalization_constant);
      new_sum += mode->weight();
    }
    printf("After renormalization sum of ratios is %g. \n", new_sum);
  }
}

}  // namespace Smash
