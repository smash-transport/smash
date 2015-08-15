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
  void reset_clock(const Clock initial_clock) {
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

/**
 * This struct contains all the additional parameters that are needed when
 * using adaptive time steps.
 */
struct AdaptiveParameters {
  float new_dt(float rate) const {
    return 0.5f * target_missed_actions / rate;
  }

  /**
   * The smoothing factor \f$ \alpha \f$ is responsible for smoothing the
   * estimate of the rate \f$ r \f$. It is used in the following equation
   * \f[ r = r_{old} + \alpha (r_{current} - r_{old}) \f]. An
   * \f$ \alpha \f$ of 1 corresponds to not taking an average at all and 0
   * corresponds to not changing the estimate of the rate at all with the new
   * data.
   */
  float smoothing_factor;

  /**
   * This is the criterion by which it is decided how long a time step should
   * be. It is the link between the estimation of the mean free time (or rate)
   * and the actual time step size that is chosen. It corresponds to the
   * number of actions per particle per time step. The value should be smaller
   * than 1 to keep the risk of missed actions small. A smaller value will
   * lead on average to smaller time steps.
   */
  float target_missed_actions;

  /**
   * This factor sets the limit by how much the criterion can be exceeded
   * before the time step is simply aborted. A value of 1 corresponds to no
   * room for overshooting at all which will lead to many aborted time steps
   * and consequently longer runtime.
   */
  float allowed_deviation;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_EXPERIMENTPARAMETERS_H_
