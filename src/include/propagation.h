/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PROPAGATION_H_
#define SRC_INCLUDE_PROPAGATION_H_

#include <vector>

class Box;
class Parameters;
class ParticleData;

void propagate_particles(std::vector<ParticleData> *particles,
  Parameters const &parameters, Box const &box);

#endif  // SRC_INCLUDE_PROPAGATION_H_
