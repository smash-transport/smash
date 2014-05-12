/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_CLOCK_H_
#define SRC_INCLUDE_CLOCK_H_

#include <cmath>

namespace Smash {

/** Clock tracks the simulation time, i.e., the time IN the simulation.
 *
 * Usage:
 * ------
 * \code
 *   Clock labtime(0.f, 0.1f);
 *   Clock endtime(10.f, 0.f);
 *   while (labtime < endtime) {
 *     // do something
 *     ++labtime;
 *   }
 * \endcode
 *
 **/
class Clock {
 public:
  /// default initializer: Timestep size is set to 0!
  Clock() : counter_(0), timestep_size_(0.f), reset_time_(0.f) {};
  /** initialize with base time and time step size.
   *
   * \param time base time
   * \param dt step size
   *
   */
  Clock(const float time, const float dt)
                : counter_(0)
                , timestep_size_(dt)
                , reset_time_(time) {}
  inline float current_time() const {
    return reset_time_ + timestep_size_ * counter_;
  }
  /// returns the time step size.
  float timestep_size() const {
    return timestep_size_;
  }
  /** sets the time step size (and resets the counter)
   *
   * \param dt new time step size
   *
   */
  void set_timestep_size(const float dt) {
    reset_time_ = current_time();
    counter_ = 0;
    timestep_size_ = dt;
  }
  /** checks if a multiple of a given interval is reached within the
   * next tick.
   *
   * \param interval The interval \f$t_i\f$ at which, for instance,
   * output is expected to happen
   *
   * \return is there a natural number n so that \f$n \cdot t_i\f$ is
   * between the current time and the next time: \f$\exists n \in
   * \mathbb{N}: t \le n \cdot t_i < t + \Delta t\f$.
   * 
   */
  bool multiple_is_in_next_tick(const float interval) const {
    // if the interval is less than or equal to the time step size, one
    // multiple will surely be within the next tick!
    if (interval <= timestep_size_) {
      return true;
    }
    // else, let's first see what n is in question:
    float n = floor((current_time() + timestep_size_) / interval);
    return (current_time() <= n * interval
         && n * interval < current_time() + timestep_size_);
  }
  /// advances the clock by one tick (\f$\Delta t\f$)
  Clock& operator++() {
    counter_++;
    return *this;
  }
  /// advances the clock by an arbitrary timestep
  Clock& operator+=(const float& big_timestep) {
    reset_time_ += big_timestep;
    return *this;
  }
  /// advances the clock by an arbitrary number of ticks.
  Clock& operator+=(const int& advance_several_timesteps) {
    counter_ += advance_several_timesteps;
    return *this;
  }
  /// compares the times between two clocks.
  bool operator<(const Clock& rhs) const {
    return current_time() < rhs.current_time();
  }
  /// compares the time of the clock against a fixed time.
  bool operator<(const float& time) const {
    return current_time() < time;
  }
  /// compares the time of the clock against a fixed time.
  bool operator>(const float& time) const {
    return current_time() > time;
  }
 private:
  /// clock tick. This is purely internal and will be reset when the
  /// timestep size is changed
  int counter_ = 0;
  /// The time step size \f$\Delta t\f$
  float timestep_size_ = 0.0f;
  /// The time of last reset (when counter_ was set to 0).
  float reset_time_ = 0.0f;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_CLOCK_H_
