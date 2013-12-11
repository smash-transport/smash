/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_SPHERE_H_
#define SRC_INCLUDE_SPHERE_H_

/* forward declarations */
class FourVector;

#include <stdint.h>
#include <time.h>
#include <cmath>

#include "../include/BoundaryConditions.h"
#include "../include/time.h"
#include "../include/Parameters.h"

class SphereBoundaryConditions : public BoundaryConditions
{
  public:
    /* default constructor with probable values */
    SphereBoundaryConditions(): radius(10.0f), timer_start(set_timer_start()) {}
    /* member funtions */
    inline timespec set_timer_start();
    /* special class funtions */
    virtual int evolve(Particles *particles, CrossSections *cross_sections);
    virtual void assign_params_specific(std::list<Parameters> *configuration);
    virtual void initial_conditions(Particles *particles);

  private:
    /* Sphere radius length */
    float radius;
    /* starting time of the simulation */
    timespec timer_start;
};



/* set the timer to the actual time in nanoseconds precision */
timespec inline SphereBoundaryConditions::set_timer_start(void) {
  timespec time;
  clock_gettime(&time);
  return time;
}

/* enforce periodic boundary conditions */
FourVector boundary_condition(FourVector position, const SphereBoundaryConditions &sphere);

#endif  // SRC_INCLUDE_SPHERE_H_
