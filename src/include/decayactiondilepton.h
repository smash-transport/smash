/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_DECAYACTIONDILEPTON_H_
#define SRC_INCLUDE_DECAYACTIONDILEPTON_H_

#include "decayaction.h"

namespace Smash {

class DecayActionDilepton : public DecayAction {
 public:
   DecayActionDilepton(const ParticleData &p, float time_of_execution);

   // do not perform any dilepton actions = leave it empty
   void perform(Particles *, size_t &) override {};

   /* generate_final_state stays the same for a first version
      maybe one can later add some output modifications and
      minor tweaks
    */
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYACTIONDILEPTON_H_
