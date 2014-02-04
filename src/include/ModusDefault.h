/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 * 
 *    GNU General Public License (GPLv3 or later)
 *
 */
#ifndef SRC_INCLUDE_MODUSDEFAULT_H_
#define SRC_INCLUDE_MODUSDEFAULT_H_

#include <stdint.h>
#include <list>

#include "../include/Particles.h"

/* forward declarations */
class Particles;
class CrossSections;
class ExperimentParameters;

/*
 * This is only a base class for actual Modus classes. The class will never be
 * used polymorphically, therefore it never needs virtual functions.
 *
 * Consider that there never are and never will be objects of type ModusDefault.
 *
 * This class is empty per default. You can add a function if you have a
 * function that is different in at least two subclasses and equal in at least
 * two subclasses.  Code that common to all goes into ExperimentImplementation.
 */
class ModusDefault {
 public:
  // never needs a virtual destructor

  // Missing functions for concrete Modus implementations:
  // void initial_conditions(Particles *particles);

  /**
   * Only needed for BoxModus. The default for all the other modi does nothing.
   */
  int sanity_check(Particles *p) { return 0; }

  void check_collision_geometry(Particles *particles,
                                CrossSections *cross_sections,
                                std::list<int> *collision_list,
                                size_t *rejection_conflict,
                                const ExperimentParameters &parameters);

  void propagate(Particles *particles, const ExperimentParameters &parameters);
};

#endif  // SRC_INCLUDE_MODUSDEFAULT_H_
