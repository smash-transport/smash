/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_BOX_H_
#define SRC_INCLUDE_BOX_H_

/* forward declarations */
class FourVector;

#include <stdint.h>
#include <time.h>
#include <cmath>

#include "../include/constants.h"

class box {
  public:
    /* default constructor with probable values */
    box(): steps_(10000), initial_condition_(1), length_(10.0),
      temperature_(0.1), energy_initial_(0),
      number_density_initial_(0), seed_(1),
      time_start_(set_timer_start()) {}
    /* member funtions */
    float inline length() const;
    void inline set_length(const float &LENGTH);
    int inline initial_condition() const;
    void inline set_initial_condition(const int &INITIAL_CONDITION);
    float inline energy_initial() const;
    void inline set_energy_initial(const float &energy);
    float inline number_density_initial() const;
    void inline set_number_density_inital(const float &number_density);
    float inline temperature() const;
    void inline set_temperature(const float &T);
    int64_t inline seed() const;
    void inline set_seed(const int64_t &RANDOMSEED);
    int inline steps() const;
    void inline set_steps(const int &STEPS);
    timespec inline time_start() const;
    timespec inline set_timer_start();

  private:
    /* number of steps */
    int steps_;
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
    /* initial seed for random generator */
    int64_t seed_;
    /* starting time of the simulation */
    timespec time_start_;
};

/* return the edge length */
float inline box::length(void) const {
  return length_;
}

/* set the edge length */
void inline box::set_length(const float &LENGTH) {
  length_ = LENGTH;
}

int inline box::steps(void) const {
  return steps_;
}

void inline box::set_steps(const int &STEPS) {
  steps_ = STEPS;
}

/* return the used initial condition */
int inline box::initial_condition(void) const {
  return initial_condition_;
}

/* set the initial condition */
void inline box::set_initial_condition(const int &INITIAL_CONDITION) {
  initial_condition_ = INITIAL_CONDITION;
}

float inline box::temperature(void) const {
  return temperature_;
}

void inline box::set_temperature(const float &T) {
  temperature_ = T;
}

float inline box::energy_initial(void) const {
  return energy_initial_;
}

void inline box::set_energy_initial(const float &energy) {
  energy_initial_ = energy;
}

float inline box::number_density_initial(void) const {
  return number_density_initial_;
}

void inline box::set_number_density_inital(const float &number_density) {
  number_density_initial_ = number_density;
}

int64_t inline box::seed(void) const {
  return seed_;
}

void inline box::set_seed(const int64_t &randomseed) {
  seed_ = randomseed;
}

timespec inline box::time_start(void) const {
  return time_start_;
}

timespec inline box::set_timer_start(void) {
  timespec time;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time);
  return time;
}

FourVector boundary_condition(FourVector position, const box &box);

#endif  // SRC_INCLUDE_BOX_H_
