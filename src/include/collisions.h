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

#include "../include/CrossSections.h"
#include "../include/Particles.h"

class BoundaryConditions;

/* populates collision list if collision applies */
void collision_criteria_geometry(Particles *particles,
  CrossSections *cross_sections,
  std::list<int> *collision_list, BoundaryConditions const &para, int id_a,
  int id_b, size_t *rejection_conflict);

/* does collisions according to collision table */
size_t collide_particles(Particles *particles, std::list<int> *collision_list,
                         size_t id_event, int *resonance_formations);

#endif  // SRC_INCLUDE_COLLISIONS_H_
