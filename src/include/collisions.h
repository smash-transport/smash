/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_COLLISIONS_H_
#define SRC_INCLUDE_COLLISIONS_H_

#include <cstdlib>
#include <list>
#include <map>
#include <vector>

class Parameters;
class ParticleData;
class ParticleType;

/* does_decay - does a resonance decay on this timestep? */
bool does_decay(std::vector<ParticleData> *particle,
  std::vector<ParticleType> *particle_type, std::map<int, int> *map_type,
  std::list<int> *collision_list, const Parameters &parameters, int id_res);

/* populates collision list if collision applies */
void collision_criteria_geometry(std::vector<ParticleData> *particle,
  std::vector<ParticleType> *type, std::map<int, int> *map_type,
  std::list<int> *collision_list, Parameters const &para, int id, int id_other);

/* does collisions according to collision table */
size_t collide_particles(std::vector<ParticleData> *particle,
  std::vector<ParticleType> *type, std::map<int, int> *map_type,
  std::list<int> *collision_list, size_t id_event);

#endif  // SRC_INCLUDE_COLLISIONS_H_
