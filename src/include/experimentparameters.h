/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_EXPERIMENTPARAMETERS_H_
#define SRC_INCLUDE_EXPERIMENTPARAMETERS_H_

#include "clock.h"

namespace Smash {

/**
 * Helper structure for Experiment.
 *
 * Experiment has one member of this struct. In essence the members of this
 * struct are members of Experiment, but combined in one structure for easier
 * function argument passing.
 */
struct ExperimentParameters {
  /// system clock (for simulation time keeping in the computational frame)
  Clock labclock;
  /// output clock to keep track of the next output time
  Clock outputclock;
  /// number of test particle
  int testparticles;
  /// width of gaussian Wigner density of particles
  float gaussian_sigma;
  /// distance at which gaussian is cut, i.e. set to zero, IN SIGMA (not fm)
  float gauss_cutoff_in_sigma;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_EXPERIMENTPARAMETERS_H_
