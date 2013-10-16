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

#include "../include/Laboratory.h"
#include "../include/time.h"

class Sphere : public Laboratory {
  public:
    /* default constructor with probable values */
    Sphere(): radius_(10.0f), time_start_(set_timer_start()) {}
    /* useful constructor with explicit values for laboratory */
    explicit Sphere(Laboratory lab): Laboratory(lab), radius_(10.0f),
      time_start_(set_timer_start()) {}
    /* member funtions */
    float inline radius() const;
    void inline set_radius(const float &RADIUS);
    timespec inline time_start() const;
    timespec inline set_timer_start();

  private:
    /* Sphere radius length */
    float radius_;
    /* starting time of the simulation */
    timespec time_start_;
};

/* return the edge length */
float inline Sphere::radius(void) const {
  return radius_;
}

/* set the edge length */
void inline Sphere::set_radius(const float &RADIUS) {
  radius_ = RADIUS;
}


/* return when the timer started */
timespec inline Sphere::time_start(void) const {
  return time_start_;
}

/* set the timer to the actual time in nanoseconds precision */
timespec inline Sphere::set_timer_start(void) {
  timespec time;
  clock_gettime(&time);
  return time;
}

/* enforce periodic boundary conditions */
FourVector boundary_condition(FourVector position, const Sphere &sphere);

#endif  // SRC_INCLUDE_SPHERE_H_
