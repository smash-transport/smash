/*
 *    Copyright (c) 2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_DECAYS_H_
#define SRC_INCLUDE_DECAYS_H_

#include <cstdlib>
#include <list>
#include <map>
#include <vector>

#include "include/particles.h"
#include "action.h"

namespace Smash {

class ModusDefault;

size_t decay_particles(Particles *particles, std::vector<ActionPtr> &decay_list,
                       size_t id_process);

/* resonance decay process */
int resonance_decay(Particles *particles, int particle_id);

/* 1->2 process kinematics */
int one_to_two(Particles *particles, int resonance_id, int type_a, int type_b);

/* 1->3 process kinematics */
int one_to_three(Particles *particles, int resonance_id,
                 int type_a, int type_b, int type_c);

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYS_H_
