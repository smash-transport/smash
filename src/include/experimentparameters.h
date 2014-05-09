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
  /// system clock (for simulation time keeping in the computational
  /// frame)
  Clock labclock;
  /// Time step size
  float timestep_size() const {
    return labclock.timestep_size();
  }
  /// returns if output should happen now
  bool need_intermediate_output() const {
    return labclock.multiple_is_in_next_tick(output_interval);
  }
  /// time interval between SMASH giving measurables
  const float output_interval;
  /// cross section of the elastic scattering
  const float cross_section;
  /// number of test particle
  int testparticles;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_EXPERIMENTPARAMETERS_H_
