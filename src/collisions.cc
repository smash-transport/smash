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
  std::list<int> *collision_list, box box, int id, int id_other) {
  double distance_squared, time_collision;

  /* distance criteria according to cross_section */
  distance_squared = particle_distance(&(*particle)[id],
    &(*particle)[id_other]);
  if (distance_squared >= box.cross_section() * fm2_mb * M_1_PI)
    return;

  /* check according timestep: positive and smaller */
  time_collision = collision_time(&(*particle)[id], &(*particle)[id_other]);
  if (time_collision < 0 || time_collision >= box.eps())
    return;

  /* check for minimal collision time */
  if ((*particle)[id].collision_time() > 0
        && time_collision > (*particle)[id].collision_time()) {
    printd("%g Not minimal particle %d <-> %d\n",
        (*particle)[id].position().x0(), id, id_other);
    return;
  }

  /* just collided with this particle */
  if ((*particle)[id].collision_time() == 0
      && id_other == (*particle)[id].collision_id()) {
    printd("%g Skipping particle %d <-> %d\n",
        (*particle)[id].position().x0(), id, id_other);
    return;
  }

  /* handle minimal collision time */
  if (unlikely((*particle)[id].collision_time() > 0)) {
    int not_id = (*particle)[id].collision_id();
    printd("Not colliding particle %d <-> %d\n", id, not_id);
    /* unset collision partner to zero time and unexisting id */
    (*particle)[not_id].set_collision(0.0, -1);
    /* remove any of those partners from the list */
    collision_list->remove(id);
    collision_list->remove(not_id);
    /* XXX: keep track of multiple possible collision partners */
  }

  /* setup collision partners */
  printd("distance squared particle %d <-> %d: %g \n", id, id_other,
    distance_squared);
  printd("collision time particle %d <-> %d: %g \n", id, id_other,
    time_collision);
  (*particle)[id].set_collision(time_collision, id_other);
  (*particle)[id_other].set_collision(time_collision, id);
  /* add to collision list */
  collision_list->push_back(id);
}

/* colliding_particle - particle interaction */
void collide_particles(std::vector<ParticleData> *particle, ParticleType *type,
  std::map<int, int> *map_type, std::list<int> *collision_list) {
  FourVector velocity_com;

  /* collide: 2 <-> 2 soft momenta exchange */
  for (std::list<int>::iterator id = collision_list->begin();
    id != collision_list->end(); ++id) {
    int id_other = (*particle)[*id].collision_id();
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
    boost_from_COM(&(*particle)[*id], &(*particle)[id_other],
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
