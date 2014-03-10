/*
 *    Copyright (c) 2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

/**
 * \file decays.h
 * Functions related to resonance decays.
 */

#ifndef SRC_INCLUDE_DECAYS_H_
#define SRC_INCLUDE_DECAYS_H_

#include <cstdlib>
#include <list>
#include <map>
#include <vector>

#include "include/particles.h"

namespace Smash {

class ModusDefault;

/**
 * Check resonances for decays in this time step.
 *
 * Go through the existing resonances and create a list
 * of those which will decay in this time step.
 */
void check_decays(Particles *particles, std::list<int> *decay_list,
                  const float timestep);

/**
 * Process the decay list.
 *
 * Go through the list of resonances decaying in this time step.
 * Boost to the resonance rest frame, do the decay,
 * boost back the product particles and remove the decayed particle
 * from the active particles data structure.
 */
size_t decay_particles(Particles *particles, std::list<int> *decay_list,
                       size_t id_process);


/// Given a resonance, select one of its decay modes and do the decay
int resonance_decay(Particles *particles, int particle_id);

/**
 * Kinematics of a 1-to-2 decay process.
 *
 * Given a resonance and the types of decay product particles,
 * sample the momenta and position of the products and add them
 * to the active particles data structure.
 */
int one_to_two(Particles *particles, int resonance_id, int type_a, int type_b);

/**
 * Kinematics of a 1-to-3 decay process.
 *
 * Given a resonance and the types of decay product particles,
 * sample the momenta and position of the products and add them
 * to the active particles data structure.
 */
int one_to_three(Particles *particles, int resonance_id,
                 int type_a, int type_b, int type_c);

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYS_H_
