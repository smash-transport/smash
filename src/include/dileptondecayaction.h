/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_DILEPTONDECAYACTION_H_
#define SRC_INCLUDE_DILEPTONDECAYACTION_H_

#include "decayaction.h"

namespace Smash {

class DileptonDecayAction : public DecayAction {
 public:
   DileptonDecayAction(const ParticleData &p, float time_of_execution);

   // do not perform any dilepton actions = leave it empty
   void perform(Particles *particles, size_t &id_process) override {};

   /* generate_final_state stays the same for a first version
      maybe one can later add some output modifications and
      minor tweaks
    */
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DILEPTONDECAYACTION_H_
