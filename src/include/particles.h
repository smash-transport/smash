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

void check_collision(ParticleData *particle,
  std::list<ParticleData> *collision_list, box box, int id, int number);

void collide_particles(ParticleData *particle,
  std::list<ParticleData> *collision_list);

#endif  // SRC_INCLUDE_PARTICLES_H_
