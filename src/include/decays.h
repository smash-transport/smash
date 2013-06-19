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
bool does_decay(ParticleData *particle, ParticleType *particle_type,
                std::list<int> *decay_list, const Parameters &parameters);

size_t decay_particles(std::map<int, ParticleData> *particle,
  std::vector<ParticleType> *type, std::map<int, int> *map_type,
  std::list<int> *decay_list, size_t id_process, size_t *largest_id);

#endif  // SRC_INCLUDE_DECAYS_H_
