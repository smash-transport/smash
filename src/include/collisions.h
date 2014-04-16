/*
 *    Copyright (c) 2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

/**
 * \file collisions.h
 * Functions related to two-particle collisions.
 */

#ifndef SRC_INCLUDE_COLLISIONS_H_
#define SRC_INCLUDE_COLLISIONS_H_

#include <list>

#include "include/crosssections.h"
#include "include/particles.h"
#include "action.h"

namespace Smash {

/**
 * Takes two particles and checks if they collide in this time step,
 * based on the geometrical interpretation of the cross section,
 * assigning an interaction for them in case of a positive check.
 *
 * Calculates the total cross section (elastic + resonance formation)
 * \f$\sigma_{tot}\f$ for the collision of the given two particles.
 *
 * Based on this, the following conditions must be fulfilled
 * for the collision to occur:
 *
 * 1. The distance between two particles must be less than their
 * interaction distance \f$\sqrt{\sigma_{tot}/\pi}\f$.
 * 2. The particles must reach their minimum distance within this time step.
 * 3. Neither of the particles is interacting with other particles
 * in this time step before the collision;
 * in other words, the collision time is a minimum for both particles.
 *
 * If all conditions are fulfilled, the specific interaction is randomly chosen
 * from the available interactions based on how large ratio of the total cross
 * section each process represents.
 *
 * The two particles are then assigned to collide with each other
 * at the determined time, with determined final state particles.
 * Both particles carry this information in their data structures,
 * while only the first particle is added to the collision list.
 *
 * \param[in] particles Particles in the simulation.
 * \param[in] cross_sections Parametrizations of (elastic) cross sections.
 * \param[in,out] collision_list List of particles assigned for collision
 * in this time step.
 * \param[in] timestep Time step size in the simulation.
 * \param[in] id_a ID of the first particle of the pair to be checked.
 * \param[in] id_b ID of the second particle of the pair to be checked.
 * \param[in,out] rejection_conflict Counter of how many times a collision
 * had to be rejected because a particle had an interaction
 * earlier within this time step.
 */
void collision_criteria_geometry(Particles *particles,
  CrossSections *cross_sections,
  std::list<int> *collision_list, const float timestep, int id_a,
  int id_b, size_t *rejection_conflict);

}  // namespace Smash

#endif  // SRC_INCLUDE_COLLISIONS_H_
