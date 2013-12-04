/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_BOX_H_
#define SRC_INCLUDE_BOX_H_

#include <stdint.h>
#include <time.h>
#include <cmath>

#include "../include/CrossSections.h"
#include "../include/Laboratory.h"
#include "../include/Particles.h"
#include "../include/time.h"

class BoxBoundaryConditions : public BoundaryConditions
{
public:
  /* default constructor with probable values */
    Box(): initial_condition_(1), length_(10.0f), temperature_(0.1f),
      energy_initial_(0.0f), number_density_initial_(0.0f),
      time_start_(set_timer_start()) {}
    /* access member funtions */
    inline int initial_condition() const;
    inline void set_initial_condition(int INITIAL_CONDITION);
    inline float length() const;
    inline void set_length(float LENGTH);
    inline float energy_initial() const;
    inline void set_energy_initial(float energy);
    inline float number_density_initial() const;
    inline void set_number_density_inital(float number_density);
    inline float temperature() const;
    inline void set_temperature(float T);
    inline timespec time_start() const;
    inline timespec set_timer_start();
    /* special class funtions */
    virtual int evolve(Particles *particles, CrossSections *cross_sections);
    virtual void process_config(char *path);
    virtual void assign_params_specific(Parameters *configuration);
    virtual void initial_conditions(Particles *particles);

  private:
    /* initial condition */
    int initial_condition_;
    /* Cube edge length */
    float length_;
    /* Temperature of the Boltzmann distribution for thermal initialization */
    float temperature_;
    /* initial total energy of the box */
    float energy_initial_;
    /* initial number density of the box */
    float number_density_initial_;
    /* starting time of the simulation */
    timespec time_start_;
};

/* return the used initial condition */
inline int Box::initial_condition(void) const {
  return initial_condition_;
}

/* set the initial condition */
inline void Box::set_initial_condition(int INITIAL_CONDITION) {
  initial_condition_ = INITIAL_CONDITION;
}

/* return the edge length */
float inline Box::length(void) const {
  return length_;
}

/* set the edge length */
void inline Box::set_length(float LENGTH) {
  length_ = LENGTH;
}

/* return the used IC temperature */
float inline Box::temperature(void) const {
  return temperature_;
}

/* set the IC temperature */
void inline Box::set_temperature(float T) {
  temperature_ = T;
}

/* return the IC total energy */
float inline Box::energy_initial(void) const {
  return energy_initial_;
}

/* set the IC total energy */
void inline Box::set_energy_initial(float energy) {
  energy_initial_ = energy;
}

/* return the IC number density */
float inline Box::number_density_initial(void) const {
  return number_density_initial_;
}

/* set the IC number density */
void inline Box::set_number_density_inital(float number_density) {
  number_density_initial_ = number_density;
}

/* return when the timer started */
timespec inline Box::time_start(void) const {
  return time_start_;
}

/* set the timer to the actual time in nanoseconds precision */
timespec inline Box::set_timer_start(void) {
  timespec time;
  clock_gettime(&time);
  return time;
}

#endif  // SRC_INCLUDE_BOX_H_
