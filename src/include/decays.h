/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_DECAYS_H_
#define SRC_INCLUDE_DECAYS_H_

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

#endif  // SRC_INCLUDE_DECAYS_H_
