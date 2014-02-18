/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_NUCLEUSMODUS_H_
#define SRC_INCLUDE_NUCLEUSMODUS_H_

#include <stdint.h>
#include <cmath>
#include <list>
#include <string>

#include "include/crosssections.h"
#include "include/modusdefault.h"
#include "include/nucleus.h"
#include "include/particles.h"
#include "include/parameters.h"

struct ExperimentParameters;

class NucleusModus : public ModusDefault {
 public:
  /* default constructor with probable values */
  NucleusModus() = default;

  /* special class funtions */
  void assign_params(
      std::list<Parameters> *configuration);  // TODO(mkretz) -> ctor
  void print_startup();  // TODO(mkretz): needs to be discoverable from an
                         // outside "printer"

  void initial_conditions(Particles *particles,
                          const ExperimentParameters &parameters);

  // in ModusDefault:
  // * sanity_check
  // * check_collision_geometry
  // * propagate

 private:
  Nucleus projectile_;
  Nucleus target_;
  // Center-of-mass energy of the collision
  float sqrts_ = 1.f;
  // initial separation of nuclei from origin
  float initial_displacement_ = 1.1;
  // steepness of woods-saxon distribution
  float steepness_ = .545;
};

#endif  // SRC_INCLUDE_NUCLEUSMODUS_H_
