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
#include <list>

#include "../include/CrossSections.h"
#include "../include/Modus.h"
#include "../include/Particles.h"
#include "../include/time.h"
#include "../include/Parameters.h"

class SphereModus : public Modus {
 public:
    /* default constructor with probable values */
    SphereModus(): radius_(10.0f), timer_start_(set_timer_start()) {}
    /* member funtions */
    /* special class funtions */
    virtual int evolve(Particles *particles, CrossSections *cross_sections);
    virtual void assign_params_specific(std::list<Parameters> *configuration);
    virtual void initial_conditions(Particles *particles);
//    virtual FourVector boundary_condition(FourVector position);
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

#endif  // SRC_INCLUDE_SPHERE_H_
