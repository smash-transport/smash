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

void DecayModes::add_mode(float ratio, int L, std::vector<PdgCode> pdg_list) {
  if (pdg_list.size() < 2) {
    throw InvalidDecay(
        "DecayModes::add_mode was instructed to add a decay mode with less "
        "than 2 particles. This is an invalid input.");
  }
  DecayBranch branch;
  branch.set_weight(ratio);
  branch.set_angular_momentum(L);
  branch.set_particles(std::move(pdg_list));
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
    for (std::vector<DecayBranch>::iterator mode
           = decay_modes_.begin(); mode != decay_modes_.end(); ++mode) {
      mode->set_weight(mode->weight() / renormalization_constant);
      new_sum += mode->weight();
    }
    printf("After renormalization sum of ratios is %g. \n", new_sum);
  }
}

}  // namespace Smash
