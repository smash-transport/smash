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

#include "include/ParticleData.h"
#include "include/ParticleType.h"
#include "include/box.h"
#include "include/constants.h"
#include "include/outputroutines.h"
#include "include/particles.h"

/* collision_criteria_geometry - check by geomatrical method if a collision
 *                               happens between particles
 */
void collision_criteria_geometry(std::vector<ParticleData> *particle,
  std::list<int> *collision_list, box box, int id_a, int id_b) {
  {
  /* distance criteria according to cross_section */
  double const distance_squared = particle_distance(&(*particle)[id_a],
    &(*particle)[id_b]);
  if (distance_squared >= box.cross_section() * fm2_mb * M_1_PI)
    return;
  printd("distance squared particle %d <-> %d: %g \n", id_a, id_b,
    distance_squared);
  }

  /* check according timestep: positive and smaller */
  double const time_collision = collision_time(&(*particle)[id_a],
    &(*particle)[id_b]);
  if (time_collision < 0 || time_collision >= box.eps())
    return;

  /* check for minimal collision time */
  if ((*particle)[id_a].collision_time() > 0
        && time_collision > (*particle)[id_a].collision_time()) {
    printd("%g Not minimal particle %d <-> %d\n",
        (*particle)[id_a].position().x0(), id_a, id_b);
    return;
  }

  /* just collided with this particle */
  if ((*particle)[id_a].collision_time() == 0
      && id_b == (*particle)[id_a].id_partner()) {
    printd("%g Skipping particle %d <-> %d\n",
        (*particle)[id_a].position().x0(), id_a, id_b);
    return;
  }

  /* handle minimal collision time */
  if (unlikely((*particle)[id_a].collision_time() > 0)) {
    int id_not = (*particle)[id_a].id_partner();
    printd("Not colliding particle %d <-> %d\n", id_a, id_not);
    /* unset collision partner to zero time and unexisting id */
    (*particle)[id_not].set_collision(0.0, -1);
    /* remove any of those partners from the list */
    collision_list->remove(id_a);
    collision_list->remove(id_not);
    /* XXX: keep track of multiple possible collision partners */
  }

  /* setup collision partners */
  printd("collision time particle %d <-> %d: %g \n", id_a, id_b,
    time_collision);
  (*particle)[id_b].set_collision(time_collision, id_b);
  (*particle)[id_a].set_collision(time_collision, id_a);
  /* add to collision list */
  collision_list->push_back(id_a);
}

/* colliding_particle - particle interaction */
void collide_particles(std::vector<ParticleData> *particle, ParticleType *type,
  std::map<int, int> *map_type, std::list<int> *collision_list) {
  FourVector velocity_com;

  /* collide: 2 <-> 2 soft momenta exchange */
  for (std::list<int>::iterator id = collision_list->begin();
    id != collision_list->end(); ++id) {
    int id_other = (*particle)[*id].id_partner();
    printd("particle types %s<->%s colliding %d<->%d time %g\n",
      type[(*map_type)[*id]].name().c_str(),
      type[(*map_type)[id_other]].name().c_str(), *id, id_other,
      (*particle)[*id].position().x0());
    write_oscar((*particle)[*id], (*particle)[id_other], type[(*map_type)[*id]],
      type[(*map_type)[id_other]], 1);
    printd("particle 1 momenta before: %g %g %g %g\n",
      (*particle)[*id].momentum().x0(), (*particle)[*id].momentum().x1(),
      (*particle)[*id].momentum().x2(), (*particle)[*id].momentum().x3());
    printd("particle 2 momenta before: %g %g %g %g\n",
      (*particle)[id_other].momentum().x0(),
      (*particle)[id_other].momentum().x1(),
      (*particle)[id_other].momentum().x2(),
      (*particle)[id_other].momentum().x3());

    /* exchange in center of momenta */
    boost_COM(&(*particle)[*id], &(*particle)[id_other], &velocity_com);
    momenta_exchange(&(*particle)[*id], &(*particle)[id_other],
      type[(*map_type)[*id]].mass(), type[(*map_type)[id_other]].mass());
    boost_back_COM(&(*particle)[*id], &(*particle)[id_other],
      &velocity_com);
    write_oscar((*particle)[*id], (*particle)[id_other], type[(*map_type)[*id]],
      type[(*map_type)[id_other]], -1);
    printd("particle 1 momenta after: %g %g %g %g\n",
      (*particle)[*id].momentum().x0(), (*particle)[*id].momentum().x1(),
      (*particle)[*id].momentum().x2(), (*particle)[*id].momentum().x3());
    printd("particle 2 momenta after: %g %g %g %g\n",
      (*particle)[id_other].momentum().x0(),
      (*particle)[id_other].momentum().x1(),
      (*particle)[id_other].momentum().x2(),
      (*particle)[id_other].momentum().x3());

    /* unset collision time for both particles + keep id */
    (*particle)[*id].set_collision_time(0.0);
    (*particle)[id_other].set_collision_time(0.0);
  }
  /* empty the collision table */
  collision_list->clear();
}
