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

class ModusDefault;

/* does_decay - does a resonance decay on this timestep? */
void check_decays(Particles *particles, std::list<int> *decay_list,
                  const float timestep);

size_t decay_particles(Particles *particles, std::list<int> *decay_list,
                       size_t id_process);

/* resonance decay process */
int resonance_decay(Particles *particles, int particle_id);

/* 1->2 process kinematics */
int one_to_two(Particles *particles, int resonance_id, int type_a, int type_b);

/* 1->3 process kinematics */
int one_to_three(Particles *particles, int resonance_id,
                 int type_a, int type_b, int type_c);

#endif  // SRC_INCLUDE_DECAYS_H_
