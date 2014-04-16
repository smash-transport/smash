/*
 *
 *    Copyright (c) 2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "include/collisions.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <vector>

#include "include/collisiondata.h"
#include "include/constants.h"
#include "include/crosssections.h"
#include "include/fourvector.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/particles.h"
#include "include/processbranch.h"
#include "include/resonances.h"
#include "include/random.h"

namespace Smash {

/* collision_criteria_geometry - check by geometrical method if a collision
 *                               happens between particles
 */
void collision_criteria_geometry(Particles *particles,
  CrossSections *cross_sections,
  std::list<int> *collision_list, const float timestep, int id_a,
  int id_b, size_t *rejection_conflict) {
  /* just collided with this particle */
  if (particles->data(id_a).id_process() >= 0
      && particles->data(id_a).id_process()
      == particles->data(id_b).id_process()) {
    printd("Skipping collided particle %d <-> %d at time %g due process %d\n",
        id_a, id_b, particles->data(id_a).position().x0(),
        particles->data(id_a).id_process());
    return;
  }

  /* Compute kinematic quantities needed for cross section calculations  */
  cross_sections->compute_kinematics(particles, id_a, id_b);
  /* Resonance production cross section */
  std::vector<ProcessBranch> resonance_xsections
    = resonance_cross_section(particles->data(id_a), particles->data(id_b),
      particles->type(id_a), particles->type(id_b), particles);

  /* Total cross section is elastic + resonance production  */
  /* (Ignore annihilation and total for now) */
  const float total_cross_section
    = cross_sections->elastic(particles, id_a, id_b)
    + resonance_xsections.at(0).weight();

  {
    /* distance criteria according to cross_section */
    const double distance_squared = particle_distance(
               particles->data_pointer(id_a), particles->data_pointer(id_b));
    if (distance_squared >= total_cross_section * fm2_mb * M_1_PI)
      return;
    printd("distance squared particle %d <-> %d: %g \n", id_a, id_b,
           distance_squared);
  }

  /* check according timestep: positive and smaller */
  const double time_collision = collision_time(particles->data(id_a),
    particles->data(id_b));
  if (time_collision < 0.0 || time_collision >= timestep)
    return;

  /* check for minimal collision time both particles */
  if ((particles->data(id_a).collision_time() > 0.0
        && time_collision > particles->data(id_a).collision_time())
      || (particles->data(id_b).collision_time() > 0
        && time_collision > particles->data(id_b).collision_time())) {
    printd("%g Not minimal particle %d <-> %d\n",
        particles->data(id_a).position().x0(), id_a, id_b);
    return;
  }

  /* handle minimal collision time of both particles */
  if (unlikely(particles->data(id_a).collision_time() > 0.0)) {
    int id_not = particles->data(id_a).id_partner();
    printd("Not colliding particle %d <-> %d\n", id_a, id_not);
    /* unset collision partner to zero time and unexisting id */
    if (particles->has_data(id_not))
      particles->data_pointer(id_not)->set_collision(-1, 0.0, -1);
    /* remove any of those partners from the list */
    if (id_a < id_not) {
      printd("Removing particle %d from collision list\n", id_a);
      collision_list->remove(id_a);
    } else {
      printd("Removing particle %d from collision list\n", id_not);
      collision_list->remove(id_not);
    }
    /* collect statistics of multiple possible collision partner */
    (*rejection_conflict)++;
  }
  if (unlikely(particles->data(id_b).collision_time() > 0.0)) {
    int id_not = particles->data(id_b).id_partner();
    printd("Not colliding particle %d <-> %d\n", id_b, id_not);
    /* unset collision partner to zero time and unexisting id */
    if (particles->has_data(id_not))
      particles->data_pointer(id_not)->set_collision(-1, 0.0, -1);
    /* remove any of those partners from the list */
    if (id_b < id_not) {
      printd("Removing particle %d from collision list\n", id_b);
      collision_list->remove(id_b);
    } else {
      printd("Removing particle %d from collision list\n", id_not);
      collision_list->remove(id_not);
    }
    /* collect statistics of multiple possible collision partner */
    (*rejection_conflict)++;
  }

  /* If resonance formation probability is high enough, do that,
   * otherwise do elastic collision
   */
  int interaction_type = 0;
  std::vector<int> final_particles;
  if (resonance_xsections.at(0).weight() > really_small) {
    double random_interaction = Random::canonical();
    float interaction_probability = 0.0;
    std::vector<ProcessBranch>::const_iterator resonances
      = resonance_xsections.begin();
    while (interaction_type == 0 && resonances != resonance_xsections.end()) {
      if (resonances->particle_list().size() > 1
          || resonances->particle_list().at(0) != 0) {
        interaction_probability += resonances->weight() / total_cross_section;
        if (random_interaction < interaction_probability) {
          interaction_type = 1;
          final_particles = resonances->particle_list();
        }
      }
      ++resonances;
    }
  }

  /* setup collision partners */
  printd("collision type %d particle %d <-> %d time: %g\n", interaction_type,
     id_a, id_b, time_collision);
  particles->data_pointer(id_a)->set_collision(interaction_type,
                                   time_collision, id_b, final_particles);
  particles->data_pointer(id_b)->set_collision(interaction_type,
                                   time_collision, id_a, final_particles);
  printd("collision type %d particle %d <-> %d time: %g\n",
         particles->data(id_a).process_type(),
         particles->data(id_a).id_partner(),
         particles->data(id_b).id_partner(),
         particles->data(id_a).collision_time());
  /* add to collision list */
  collision_list->push_back(id_a);
  resonance_xsections.clear();
}


}  // namespace Smash
