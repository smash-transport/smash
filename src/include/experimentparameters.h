/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_EXPERIMENTPARAMETERS_H_
#define SRC_INCLUDE_EXPERIMENTPARAMETERS_H_

#include "include/clock.h"

namespace Smash {

/**
 * Helper structure for Experiment.
 *
 * Experiment has one member of this struct. In essence the members of this
 * struct are members of Experiment, but combined in one structure for easier
 * function argument passing.
 */
struct ExperimentParameters {
  /// temporal time step
  Clock labclock;
  /// cross section of the elastic scattering
  const float cross_section;
  /// number of test particle
  int testparticles;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_EXPERIMENTPARAMETERS_H_
