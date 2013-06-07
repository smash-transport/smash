/*
 *
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include "include/collisions.h"

#include <cstdio>

#include "include/Parameters.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"
#include "include/constants.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/particles.h"

/* collision_criteria_geometry - check by geometrical method if a collision
 *                               happens between particles
 */
void collision_criteria_geometry(std::vector<ParticleData> *particle,
  std::vector<ParticleType> *particle_type, std::map<int, int> *map_type,
  std::list<int> *collision_list, const Parameters &parameters, int id_a,
  int id_b) {

  /* just collided with this particle */
  if ((*particle)[id_a].id_process() >= 0
      && (*particle)[id_a].id_process() == (*particle)[id_b].id_process()) {
    printd("Skipping collided particle %d <-> %d at time %g\n",
        id_a, id_b, (*particle)[id_a].position().x0());
    return;
  }

  /* Resonance production cross section */
  const double resonance_xsection = resonance_cross_section(
   &(*particle)[id_a], &(*particle)[id_b], &(*particle_type)[(*map_type)[id_a]],
   &(*particle_type)[(*map_type)[id_b]], particle_type);

  /* Total cross section is elastic + resonance production  */
  const double total_cross_section = parameters.cross_section()
    + resonance_xsection;

  {
    /* distance criteria according to cross_section */
    const double distance_squared = particle_distance(&(*particle)[id_a],
                                     &(*particle)[id_b]);
    if (distance_squared >= total_cross_section * fm2_mb * M_1_PI)
      return;
    printd("distance squared particle %d <-> %d: %g \n", id_a, id_b,
           distance_squared);
  }

  /* check according timestep: positive and smaller */
  const double time_collision = collision_time(&(*particle)[id_a],
    &(*particle)[id_b]);
  if (time_collision < 0 || time_collision >= parameters.eps())
    return;

  /* check for minimal collision time both particles */
  if (((*particle)[id_a].collision_time() > 0
        && time_collision > (*particle)[id_a].collision_time())
      || ((*particle)[id_b].collision_time() > 0
        && time_collision > (*particle)[id_b].collision_time())) {
    printd("%g Not minimal particle %d <-> %d\n",
        (*particle)[id_a].position().x0(), id_a, id_b);
    return;
  }

  /* handle minimal collision time of both particles */
  /* XXX: keep track of multiple possible collision partners */
  if (unlikely((*particle)[id_a].collision_time() > 0)) {
    int id_not = (*particle)[id_a].id_partner();
    printd("Not colliding particle %d <-> %d\n", id_a, id_not);
    /* unset collision partner to zero time and unexisting id */
    (*particle)[id_not].set_collision(-1, 0.0, -1);
    /* remove any of those partners from the list */
    if (id_a < id_not) {
      printd("Removing particle %d from collision list\n", id_a);
      collision_list->remove(id_a);
    } else {
      printd("Removing particle %d from collision list\n", id_not);
      collision_list->remove(id_not);
    }
  }
  if (unlikely((*particle)[id_b].collision_time() > 0)) {
    int id_not = (*particle)[id_b].id_partner();
    printd("Not colliding particle %d <-> %d\n", id_b, id_not);
    /* unset collision partner to zero time and unexisting id */
    (*particle)[id_not].set_collision(-1, 0.0, -1);
    /* remove any of those partners from the list */
    if (id_b < id_not) {
      printd("Removing particle %d from collision list\n", id_b);
      collision_list->remove(id_b);
    } else {
      printd("Removing particle %d from collision list\n", id_not);
      collision_list->remove(id_not);
    }
  }

  /* If resonance formation probability is high enough, do that,
   * otherwise do elastic collision
   */
  int interaction_type = 0;
  if (drand48() < resonance_xsection / total_cross_section)
    interaction_type = 1;

  /* setup collision partners */
  printd("collision time particle %d <-> %d: %g \n", id_a, id_b,
    time_collision);
  (*particle)[id_a].set_collision(interaction_type, time_collision, id_b);
  (*particle)[id_b].set_collision(interaction_type, time_collision, id_a);
  /* add to collision list */
  collision_list->push_back(id_a);
}

/* colliding_particle - particle interaction */
size_t collide_particles(std::vector<ParticleData> *particle,
  std::vector<ParticleType> *type, std::map<int, int> *map_type,
  std::list<int> *collision_list, size_t id_process) {
  FourVector velocity_CM;

  /* XXX: print debug output of collision list */
  /* collide: 2 <-> 2 soft momenta exchange */
  for (std::list<int>::iterator id = collision_list->begin();
    id != collision_list->end(); ++id) {
    /* relevant particle id's for the collision */
    int id_a = *id;
    int id_b = 0;
    int interaction_type = (*particle)[id_a].process_type();
    if (interaction_type < 2) {
      id_b = (*particle)[id_a].id_partner();
      printd("Process %lu type %i particle %s<->%s colliding %d<->%d time %g\n",
        id_process, interaction_type, (*type)[(*map_type)[id_a]].name().c_str(),
             (*type)[(*map_type)[id_b]].name().c_str(), id_a, id_b,
             (*particle)[id_a].position().x0());
      write_oscar((*particle)[id_a], (*particle)[id_b],
                  (*type)[(*map_type)[id_a]], (*type)[(*map_type)[id_b]], 1);
      printd("particle 1 momenta before: %g %g %g %g\n",
          (*particle)[id_a].momentum().x0(), (*particle)[id_a].momentum().x1(),
          (*particle)[id_a].momentum().x2(), (*particle)[id_a].momentum().x3());
      printd("particle 2 momenta before: %g %g %g %g\n",
          (*particle)[id_b].momentum().x0(), (*particle)[id_b].momentum().x1(),
          (*particle)[id_b].momentum().x2(), (*particle)[id_b].momentum().x3());

      /* processes computed in the center of momenta */
      boost_CM(&(*particle)[id_a], &(*particle)[id_b], &velocity_CM);
    }

    if (interaction_type == 0) {
      /* 2->2 elastic scattering*/
      printd("Process: Elastic collision.\n");
      momenta_exchange(&(*particle)[id_a], &(*particle)[id_b],
       (*type)[(*map_type)[id_a]].mass(), (*type)[(*map_type)[id_b]].mass());
    } else if (interaction_type == 1) {
      /* 2->1 resonance formation */
      printd("Process: Resonance formation.\n");
      resonance_formation(particle, type, map_type, &id_a, &id_b);
    } else if (interaction_type == 2) {
      /* 1->2 resonance decay */
      printd("Process: Resonance decay.\n");
      /* XXX: Need to boost to rest frame */
      resonance_decay(particle, type, map_type, &id_a);
      /* XXX: Boost the new particles to computational frame */
    } else
       printf("Warning: Unspecified process type, nothing done.\n");

    if (interaction_type < 2) {
      boost_back_CM(&(*particle)[id_a], &(*particle)[id_b],
       &velocity_CM);
      write_oscar((*particle)[id_a], (*particle)[id_b],
       (*type)[(*map_type)[id_a]], (*type)[(*map_type)[id_b]], -1);
      printd("particle 1 momenta after: %g %g %g %g\n",
       (*particle)[id_a].momentum().x0(), (*particle)[id_a].momentum().x1(),
       (*particle)[id_a].momentum().x2(), (*particle)[id_a].momentum().x3());
      printd("particle 2 momenta after: %g %g %g %g\n",
       (*particle)[id_b].momentum().x0(), (*particle)[id_b].momentum().x1(),
       (*particle)[id_b].momentum().x2(), (*particle)[id_b].momentum().x3());
      /* unset collision time for both particles + keep id + unset partner */
      (*particle)[id_b].set_collision_past(id_process);
    }
    (*particle)[id_a].set_collision_past(id_process);
    id_process++;
  }
  /* empty the collision table */
  collision_list->clear();
  /* return how many processes we handled */
  return id_process;
}
