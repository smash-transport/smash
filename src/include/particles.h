/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PARTICLES_H_
#define SRC_INCLUDE_PARTICLES_H_

#include <list>
class ParticleData;
class box;

/* populates collision list if collision applies */
void check_collision_criteria(ParticleData *particle,
  std::list<int> *collision_list, box box, int id, int id_other);

/* does collisions according to collision table */
void collide_particles(ParticleData *particle,
  std::list<int> *collision_list);

#endif  // SRC_INCLUDE_PARTICLES_H_
