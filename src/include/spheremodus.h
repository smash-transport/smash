/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SPHEREMODUS_H_
#define SRC_INCLUDE_SPHEREMODUS_H_

#include <stdint.h>
#include <time.h>
#include <cmath>
#include <list>

#include "include/crosssections.h"
#include "include/modusdefault.h"
#include "include/particles.h"
#include "include/time.h"

namespace Smash {

/* forward declarations */
class FourVector;

class SphereModus : public ModusDefault {
 public:
  /* default constructor with probable values */
  explicit SphereModus(Configuration modus_config);
  /* member funtions */
  /* special class funtions */
  int evolve(Particles *particles, CrossSections *cross_sections);
  float initial_conditions(Particles *particles);
  //     FourVector boundary_condition(FourVector position);

 private:
  /* Sphere radius length */
  float radius_;
};

/* enforce periodic boundary conditions */
// FourVector SphereModus::boundary_condition(FourVector position);

}  // namespace Smash

#endif  // SRC_INCLUDE_SPHEREMODUS_H_
