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

#include "include/configuration.h"
#include "include/crosssections.h"
#include "include/modusdefault.h"
#include "include/nucleus.h"
#include "include/particles.h"
#include "include/parameters.h"

struct ExperimentParameters;

class NucleusModus : public ModusDefault {
 public:
  /* default constructor with probable values */
  NucleusModus(Configuration& config);

  void initial_conditions(Particles *particles,
                          const ExperimentParameters &parameters);

 private:
  Nucleus projectile_;
  Nucleus target_;
  /// Center-of-mass energy of the individual nucleon-nucleon collisions.
  /// Note that sqrt(s_pp) != sqrt(s_nn) [i.e., sqrt(s) depends on the
  /// masses of the pair of particles we are looking at), therefore, we
  /// also need to specify which particle pair this really refers to.
  /// That is done with pdg_sNN_1_ and pdg_sNN_2_.
  float sqrt_s_NN_ = 1.f;
  /// PDG code of the first particle in the pair that defines sqrt(s_NN)
  int pdg_sNN_1_ = 2212;
  /// PDG code of the second particle in the pair that defines sqrt(s_NN)
  int pdg_sNN_2_ = 2212;
};

#endif  // SRC_INCLUDE_NUCLEUSMODUS_H_
