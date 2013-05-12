/*
 *    Copyright (c) 2012 maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_INITIAL_CONDITIONS_H_
#define SRC_INCLUDE_INITIAL_CONDITIONS_H_

#include <map>

/* forward declarations */
class ParticleData;
class ParticleType;
class box;

/* initialisation functions */
ParticleType* initial_particles(ParticleType *type);
ParticleData* initial_conditions(ParticleData *particles,
  ParticleType *particle_type, std::map<int, int> *map_type, int &number,
  box *box);

#endif  // SRC_INCLUDE_INITIAL_CONDITIONS_H_
