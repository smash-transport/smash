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
                                         float time_of_execution)
    : DecayAction({p}, time_of_execution) {}

}  // namespace Smash
