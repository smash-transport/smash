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
  /// system clock (for simulation time keeping in the computational
  /// frame)
  Clock labclock;
  /// Make one timestep to reach given time
  void reach_time_in_one_timestep(float time) {
    labclock.set_timestep_duration(time - labclock.current_time());
    ++labclock;
  }
  /// Time step size
  float timestep_duration() const {
    return labclock.timestep_duration();
  }
  /// returns if output should happen now
  bool need_intermediate_output() const {
    return labclock.multiple_is_in_next_tick(output_interval);
  }
  /// sets the time step such that it ends on the next output time
  void set_timestep_for_next_output() {
    labclock.end_tick_on_multiple(output_interval);
  }
  /// returns if the current time is exactly an output time
  bool is_output_time() const {
    return need_intermediate_output() &&
           labclock.next_time() < labclock.next_multiple(output_interval);
  }
  /// replaces the current clock with a new one.
  void reset_clock(const Clock &initial_clock) {
    labclock = std::move(initial_clock);
  }
  /// this is the time particles will have after propagating through the
  /// current time step.
  float new_particle_time() const {
    // I'm not certain if they should have the current time or the next
    // tick.
    return labclock.current_time();
        // labclock.next_time();
  }
  /// time interval between SMASH giving measurables
  const float output_interval;
  /// number of test particle
  int testparticles;
  /// width of gaussian Wigner density of particles
  float gaussian_sigma;
  /// distance at which gaussian is cut, i.e. set to zero, IN SIGMA (not fm)
  float gauss_cutoff_in_sigma;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_EXPERIMENTPARAMETERS_H_
