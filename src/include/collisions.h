/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_COLLISIONS_H_
#define SRC_INCLUDE_COLLISIONS_H_

#include <list>
#include <map>
#include <vector>

class ParticleData;
class ParticleType;
class box;

/* populates collision list if collision applies */
void collision_criteria_geometry(std::vector<ParticleData> *particle,
  std::list<int> *collision_list, box box, int id, int id_other);

/* does collisions according to collision table */
void collide_particles(std::vector<ParticleData> *particle, ParticleType *type,
  std::map<int, int> *map_type, std::list<int> *collision_list);

#endif  // SRC_INCLUDE_COLLISIONS_H_
