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

namespace Smash {

/** Clock tracks the simulation time, i.e., the time IN the simulation.
 *
 * Usage:
 * ------
 * \code
 * ...
 * \endcode
 *
 **/

class Clock {
 public:
  Clock();
  Clock(const float time, const float dt)
                : counter_(0)
                , timestep_size_(dt)
                , reset_time_(time) {}
  inline float current_time() const {
    return reset_time_ + timestep_size_ * counter_;
  }
  float timestep_size() const {
    return timestep_size_;
  }
  void set_timestep_size(const float dt) {
    reset_time_ = current_time();
    counter_ = 0;
    timestep_size_ = dt;
  }
  Clock& operator++() {
    counter_++;
    return *this;
  }
  Clock& operator+=(const float& big_timestep) {
    reset_time_ += big_timestep;
    return *this;
  }
  Clock& operator+=(const int& advance_several_timesteps) {
    counter_ += advance_several_timesteps;
    return *this;
  }
  bool operator<(const Clock& rhs) {
    return current_time() < rhs.current_time();
  }
 private:
  int counter_ = 0; // no accessor
  float timestep_size_ = 0.0; // accessor
  float reset_time_ = 0.0; // no accessor
};

}  // namespace Smash

#endif  // SRC_INCLUDE_CLOCK_H_
