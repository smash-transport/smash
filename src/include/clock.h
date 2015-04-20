/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_CLOCK_H_
#define SRC_INCLUDE_CLOCK_H_

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>

#include "logging.h"

namespace Smash {

/** Clock tracks the simulation time, i.e., the time IN the simulation.
 *
 * The basic unit is 1 fm/c = \f$1 / 2.99798542 \cdot 10^{-23}\f$s
 * \f$\approx 0.33 \cdot 10^{-24}\f$ s.
 * The resolution of the clock is 0.001 fm/c. I.e. only multiples of 0.001 fm/c
 * are representable internally.
 *
 * Usage:
 * ------
 * \code
 *   Clock labtime(0.f, 0.1f);
 *   Clock endtime(10.f, 0.f);
 *   while (labtime < endtime) {
 *     // do something
 *     // adapt the timestep size to external circumstances:
 *     if (system_is_very_dense()) {
 *       labtime.set_timestep_duration(labtime.timestep_duration() / 2.0f);
 *     }
 *     if (system_is_very_dilute()) {
 *       labtime.set_timestep_duration(labtime.timestep_duration() * 2.0f);
 *     }
 *     // let the clock tick
 *     ++labtime;
 *   }
 * \endcode
 *
 * Possible actions for Clock are:
 * \li look at it and find out the current time
 * \see current_time()
 * \see next_time()
 * \li advance the clock (by one tick, by several ticks, or by a given
 * time)
 * \see operator++()
 * \see operator+=(T)
 * \see operator+=(Representation)
 * \li set / retrieve the timestep (length of one tick)
 * \see set_timestep_duration() \see timestep_duration()
 * \li compare time against different clock or fixed value
 * \see operator<(const Clock&) const
 * \see operator<(float) const
 * \see operator>(float) const
 *
 * Internals
 * ---------
 *
 * Clock stores a time step size \f$\Delta t\f$ and a base time
 * \f$t_0\f$ as well as a counter \f$n\f$. The current time is
 * calculated from \f$t = t_0 + n \cdot \Delta t\f$. When \f$\Delta t\f$
 * is changed, \f$t_0\f$ is reset.
 *
 **/
class Clock {
  /**
   * Defines the resolution of the clock (i.e. the smallest representable time
   * difference).
   *
   * \fpPrecision
   * This constant is only used to derive other constant expressions, which are
   * in single-precision:
   * \li It doesn't matter for the execution / data-representation of SMASH at
   * all
   * \li The increased precision is necessery to determine the correctly rounded
   * inverse for \c from_float.
   */
  static constexpr double resolution = 0.00001;
  static constexpr float to_float = static_cast<float>(resolution);
  static constexpr float from_float = static_cast<float>(1. / resolution);

 public:
  /// The type used for counting ticks/time.
  using Representation = std::int64_t;
  /// default initializer: Timestep size is set to 0!
  Clock() = default;
  /** initialize with base time and time step size.
   *
   * \param time base time
   * \param dt step size
   *
   */
  Clock(const float time, const float dt)
      : timestep_duration_(dt * from_float), reset_time_(time * from_float) {
    if (dt < 0.f) {
      throw std::range_error("No negative time increment allowed");
    }
  }
  /// returns the current time
  float current_time() const {
    return (reset_time_ + timestep_duration_ * counter_) * to_float;
  }
  /** returns the time in the next tick
   *
   * This function is needed, because current_time() + timestep_duration()
   * is not the same as the next tick (numerically; this is due to
   * floating point arithmetic).
   */
  float next_time() const {
    if (counter_ * timestep_duration_ >=
        std::numeric_limits<Representation>::max() - timestep_duration_) {
      throw std::overflow_error("Too many timesteps, clock overflow imminent");
    }
    return (reset_time_ + timestep_duration_ * (counter_ + 1)) * to_float;
  }
  /// returns the time step size.
  float timestep_duration() const { return timestep_duration_ * to_float; }
  /** sets the time step size (and resets the counter)
   *
   * \param dt new time step size
   *
   */
  void set_timestep_duration(const float dt) {
    if (dt < 0.f) {
      throw std::range_error("No negative time increment allowed");
    }
    reset_time_ += timestep_duration_ * counter_;
    counter_ = 0;
    timestep_duration_ = dt * from_float;
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
   * Internally, it checks if \f$n\f$ is the same for this time step as
   * for the next time step. If so, the next multiple is outside the
   * current time step, if not, then the multiple must be within.
   *
   */
  bool multiple_is_in_next_tick(const float interval) const {
    if (interval < 0.f) {
      throw std::range_error("Negative interval makes no sense for clock");
    }
    const Representation int_interval = interval * from_float;
    // if the interval is less than or equal to the time step size, one
    // multiple will surely be within the next tick!
    if (int_interval <= timestep_duration_) {
      return true;
    }
    const auto next = reset_time_ + timestep_duration_ * counter_;
    if (unlikely(next < 0)) {
      return -next % int_interval < timestep_duration_;
    }
    return (timestep_duration_ - 1 + next) % int_interval < timestep_duration_;
  }
  /** returns the next multiple of a given interval
   *
   * \param interval the interval in question
   *
   * \return The smallest multiple of \p interval that is larger than
   * the current time.
   */
  float next_multiple(const float interval) const {
    const Representation int_interval = interval * from_float;
    const auto current = reset_time_ + timestep_duration_ * counter_;
    if (unlikely(current < 0)) {
      return current / int_interval * int_interval * to_float;
    }
    return (current / int_interval + 1) * int_interval * to_float;
  }
  /** resets the time to a pre-defined value
   *
   * This is the only way of turning the clock back. It is needed so
   * that the time can be adjusted after initialization (different
   * initial conditions may require different starting times).
   *
   **/
  void reset(const float reset_time) {
    if (reset_time < current_time()) {
      logger<LogArea::Clock>().info("Resetting clock from", current_time(),
                                    " fm/c to ", reset_time, " fm/c");
    }
    reset_time_ = reset_time * from_float;
    counter_ = 0;
  }
  /** advances the clock by one tick (\f$\Delta t\f$)
   *
   * This operator is used as `++clock`. The operator `clock++` is not
   * implemented deliberately, because that requires a copy of the clock
   * being created.
   */
  Clock& operator++() {
    // guard against overflow:
    if (counter_ * timestep_duration_ >=
        std::numeric_limits<Representation>::max() - timestep_duration_) {
      throw std::overflow_error("Too many timesteps, clock overflow imminent");
    }
    ++counter_;
    return *this;
  }
  /**
   * Advances the clock by an arbitrary timestep (multiple of 0.001 fm/c).
   *
   * \note It uses a template parameter only for disambiguation with the
   * overload below.
   */
  template <typename T>
  typename std::enable_if<std::is_floating_point<T>::value, Clock&>::type
  operator+=(T big_timestep) {
    if (big_timestep < 0.f) {
      throw std::range_error("The clock cannot be turned back.");
    }
    reset_time_ += big_timestep * from_float;
    return *this;
  }
  /// advances the clock by an arbitrary number of ticks.
  Clock& operator+=(Representation advance_several_timesteps) {
    if (counter_ * timestep_duration_ >=
        std::numeric_limits<Representation>::max() -
            advance_several_timesteps * timestep_duration_) {
      throw std::overflow_error("Too many timesteps, clock overflow imminent");
    }
    counter_ += advance_several_timesteps;
    return *this;
  }
  /// compares the times between two clocks.
  bool operator<(const Clock& rhs) const {
    return (reset_time_ + timestep_duration_ * counter_) <
           (rhs.reset_time_ + rhs.timestep_duration_ * rhs.counter_);
  }
  /// compares the time of the clock against a fixed time.
  bool operator<(float time) const { return current_time() < time; }
  /// compares the time of the clock against a fixed time.
  bool operator>(float time) const { return current_time() > time; }

 private:
  /// clock tick. This is purely internal and will be reset when the
  /// timestep size is changed
  Representation counter_ = 0;
  /// The time step size \f$\Delta t\f$ in $10^{-3}$ fm.
  Representation timestep_duration_ = 0u;
  /// The time of last reset (when counter_ was set to 0).
  Representation reset_time_ = 0;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_CLOCK_H_
