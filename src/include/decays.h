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

#include "../include/Particles.h"

class Parameters;

/* does_decay - does a resonance decay on this timestep? */
void check_decays(Particles *particles, std::list<int> *decay_list,
                  const Parameters &parameters);

size_t decay_particles(Particles *particles, std::list<int> *decay_list,
                       size_t id_process);

#endif  // SRC_INCLUDE_DECAYS_H_
