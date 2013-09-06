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

#include "../include/Particles.h"

class Parameters;

/* populates collision list if collision applies */
void collision_criteria_geometry(Particles *particles,
  std::list<int> *collision_list, Parameters const &para, int id_a,
  int id_b, size_t *rejection_conflict);

/* does collisions according to collision table */
size_t collide_particles(Particles *particles, std::list<int> *collision_list,
                         size_t id_event, int *resonance_formations);

#endif  // SRC_INCLUDE_COLLISIONS_H_
