/*
 *    Copyright (c) 2012 maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_INITIAL_CONDITIONS_H_
#define SRC_INCLUDE_INITIAL_CONDITIONS_H_

#include <map>
#include <vector>

/* forward declarations */
class Box;
class Parameters;
class ParticleData;
class ParticleType;

/* initialisation functions */
ParticleType* initial_particles(ParticleType *type);
void initial_conditions(std::vector<ParticleData> *particles,
  ParticleType *particle_type, std::map<int, int> *map_type,
  Parameters *parameters, Box *box);

#endif  // SRC_INCLUDE_INITIAL_CONDITIONS_H_
