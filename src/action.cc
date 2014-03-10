/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/action.h"

namespace Smash {

Action::Action(const ParticleList &ingoing_particles_ids,
               float time_of_execution)
    : ingoing_particles_ids_(ingoing_particles_ids),
      time_of_execution_(time_of_execution) {}

Action::~Action() {}


}  // namespace Smash
