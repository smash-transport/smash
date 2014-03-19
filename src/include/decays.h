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
 *
 * \param[in] particles Particles in the simulation.
 * \param[out] decay_list List of particles assigned for decay
 * in this time step.
 * \param[in] timestep
 *
 */
void check_decays(Particles *particles, std::list<int> *decay_list,
                  const float timestep);

/**
 * Process the decay list.
 *
 * Go through the list of resonances decaying in this time step.
 * Boost to the resonance rest frame, do the decay
 * (by calling resonance_decay),
 * boost back the product particles and remove the decayed particle
 * from the active particles data structure.
 *
 * \param[in,out] particles Particles in the simulation.
 * \param[in] decay_list List of particles assigned for decay
 * in this time step.
 * \param[in] id_process Process ID of the latest interaction
 * in the simulation; indicates the number of interactions processed
 * so far.
 *
 * \return The number of processes recorded so far in the simulation.
 */
size_t decay_particles(Particles *particles, std::list<int> *decay_list,
                       size_t id_process);

/**
 * Execute a decay process for the selected particle.
 *
 * Randomly select one of the decay modes of the particle
 * according to their relative weights. Then decay the particle
 * by calling function one_to_two or one_to_three.
 *
 * \param[in,out] particles Particles in the simulation.
 * \param[in] particle_id ID of the particle which should decay.
 *
 * \return The ID of the first decay product.
 */
int resonance_decay(Particles *particles, int particle_id);

/**
 * Kinematics of a 1-to-2 decay process.
 *
 * Given a resonance ID and the PDG codes of decay product particles,
 * sample the momenta and position of the products and add them
 * to the active particles data structure.
 *
 * \param[in,out] particles Particles in the simulation.
 * \param[in] resonance_id ID of the resonance.
 * \param[in] type_a PDG code of the first decay product.
 * \param[in] type_b PDG code of the second decay product.
 *
 * \return The ID of the first new particle.
 */
int one_to_two(Particles *particles, int resonance_id, int type_a, int type_b);

/**
 * Kinematics of a 1-to-3 decay process.
 *
 * Given a resonance ID and the PDG codes of decay product particles,
 * sample the momenta and position of the products and add them
 * to the active particles data structure.
 *
 * \param[in,out] particles Particles in the simulation.
 * \param[in] resonance_id ID of the resonance.
 * \param[in] type_a PDG code of the first decay product.
 * \param[in] type_b PDG code of the second decay product.
 * \param[in] type_c PDG code of the third decay product.
 *
 * \return The ID of the first new particle.
 */
int one_to_three(Particles *particles, int resonance_id,
                 int type_a, int type_b, int type_c);

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYS_H_
