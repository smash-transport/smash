/*
 *    Copyright (c) 2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SPHEREMODUS_H_
#define SRC_INCLUDE_SPHEREMODUS_H_

/* forward declarations */
class FourVector;

#include <stdint.h>
#include <time.h>
#include <cmath>
#include <list>

#include "include/crosssections.h"
#include "include/modusdefault.h"
#include "include/particles.h"
#include "include/time.h"
#include "include/parameters.h"

class SphereModus : public ModusDefault {
 public:
  /* default constructor with probable values */
  SphereModus(Configuration &config);
  /* member funtions */
  /* special class funtions */
  int evolve(Particles *particles, CrossSections *cross_sections);
  void initial_conditions(Particles *particles);
  //     FourVector boundary_condition(FourVector position);
  inline timespec set_timer_start();

 private:
  /* Sphere radius length */
  float radius_;
  /* starting time of the simulation */
  timespec timer_start_;
};

/* set the timer to the actual time in nanoseconds precision */
timespec inline SphereModus::set_timer_start(void) {
  timespec time;
  clock_gettime(&time);
  return time;
}

/* enforce periodic boundary conditions */
// FourVector SphereModus::boundary_condition(FourVector position);

#endif  // SRC_INCLUDE_SPHEREMODUS_H_
