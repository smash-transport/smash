/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_BOX_H_
#define SRC_INCLUDE_BOX_H_

#include <stdint.h>

class box {
  public:
    /* default constructor with improbable values */
    box(): steps_(0) {}
    /* member funtions */
    void inline set(const float A, const float EPS, int64_t RANDOMSEED,
      const float sigma, const int STEPS, const float T, const int UPDATE);
    float inline a() const;
    void inline set_a(const float &A);
    float inline cross_section() const;
    void inline set_cross_section(const float &sigma);
    float inline energy_initial() const;
    void inline set_energy_initial(const float &energy);
    float inline eps() const;
    void inline set_eps(const float &EPS);
    float inline temperature() const;
    void inline set_temperature(const float &T);
    int64_t inline seed() const;
    void inline set_seed(const int64_t &RANDOMSEED);
    int inline steps() const;
    void inline set_steps(const int &STEPS);
    int inline update() const;
    void inline set_update(const int &UPDATE);

  private:
    /* number of steps */
    int steps_;
    /* number of steps before giving measurables */
    int update_;
    /* Cube edge length */
    float a_;
    /* temporal time step */
    float eps_;
    /* Temperature of the Boltzmann distribution for thermal initialization */
    float temperature_;
    /* cross section of the elastic scattering */
    float cross_section_;
    /* initial total energy of the box */
    float energy_initial_;
    /* initial seed for random generator */
    int64_t seed_;
};

void inline box::set(const float A, const float EPS, const int64_t RANDOMSEED,
      const float sigma, const int STEPS, const float T, const int UPDATE) {
  a_ = A;
  cross_section_ = sigma;
  eps_ = EPS;
  seed_ = RANDOMSEED;
  steps_ = STEPS;
  update_ = UPDATE;
  temperature_ = T;
}

float inline box::a(void) const {
  return a_;
}

void inline box::set_a(const float &A) {
  a_ = A;
}

float inline box::eps(void) const {
  return eps_;
}

void inline box::set_eps(const float &EPS) {
  eps_ = EPS;
}

int inline box::steps(void) const {
  return steps_;
}

void inline box::set_steps(const int &STEPS) {
  steps_ = STEPS;
}

int inline box::update(void) const {
  return update_;
}

void inline box::set_update(const int &UPDATE) {
  update_ = UPDATE;
}

float inline box::temperature(void) const {
  return temperature_;
}

void inline box::set_temperature(const float &T) {
  temperature_ = T;
}

float inline box::cross_section(void) const {
  return cross_section_;
}

void inline box::set_cross_section(const float &sigma) {
  cross_section_ = sigma;
}

float inline box::energy_initial(void) const {
  return energy_initial_;
}

void inline box::set_energy_initial(const float &energy) {
  energy_initial_ = energy;
}

int64_t inline box::seed(void) const {
  return seed_;
}

void inline box::set_seed(const int64_t &randomseed) {
  seed_ = randomseed;
}

/* support for gcc branch prediction */
#ifdef __GNUC__
#define likely(x)       __builtin_expect((x), 1)
#define unlikely(x)     __builtin_expect((x), 0)
#else
#define likely(x)       (x)
#define unlikely(x)     (x)
#endif

#endif  // SRC_INCLUDE_BOX_H_
