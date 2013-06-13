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
void initial_particles(std::vector<ParticleType> *type);
void initial_conditions(std::map<int, ParticleData> *particles,
  std::vector <ParticleType> *particle_type, std::map<int, int> *map_type,
  Parameters *parameters, Box *box);

#endif  // SRC_INCLUDE_INITIAL_CONDITIONS_H_
