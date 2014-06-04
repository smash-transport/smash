/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/particledata.h"

#include "include/particletype.h"
#include "include/processbranch.h"
#include "include/decaymodes.h"
#include "include/width.h"

namespace Smash {

float ParticleData::total_width() const {
  float w = 0., partial_width_at_pole;
  if (type_->is_stable()) return w;
  const std::vector<DecayBranch> decaymodes
        = DecayModes::find(pdgcode()).decay_mode_list();
  // loop over decay modes and sum up all partial widths
  for (std::vector<DecayBranch>::const_iterator mode = decaymodes.begin();
       mode != decaymodes.end(); ++mode) {
    partial_width_at_pole = type_->width_at_pole()*mode->weight();
    if (mode->pdg_list().size()==2) {
      // mass-dependent width for 2-body decays
      w = w + width_Manley (momentum_.abs(), type_->mass(),
                            ParticleType::find(mode->pdg_list()[0]).mass(),
                            ParticleType::find(mode->pdg_list()[1]).mass(),
                            mode->angular_momentum(),partial_width_at_pole);
    }
    else {
      // constant width for three-body decays
      w = w + partial_width_at_pole;
    }
  }

  return w;
}


}  // namespace Smash
