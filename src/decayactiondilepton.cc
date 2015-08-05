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
                                         float time,
                                         float shining_weight)
    : DecayAction({p}, time), shining_weight_(shining_weight) {}

}  // namespace Smash
