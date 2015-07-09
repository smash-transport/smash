/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */


#include "include/decayactiondilepton.h"

namespace Smash {

DecayActionDilepton::DecayActionDilepton(const ParticleData &p,
                                         float time_of_execution,
                                         float shining_weight)
    : DecayAction({p}, time_of_execution), shining_weight_(shining_weight) {}


void DecayActionDilepton::dalitz_cinematics() {

}

void DecayActionDilepton::generate_final_state() {

}


}  // namespace Smash
