/*
 *
 *    Copyright (c) 2014-2018
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

namespace smash {

/**
 * Clock tracks the time in the simulation.
 *
 * The basic unit is 1 fm/c = \f$1 / 2.99798542 \cdot 10^{-23}\f$s
 * \f$\approx 0.33 \cdot 10^{-24}\f$ s.
 * The resolution of the clock is 0.000001 fm/c. I.e. only multiples of 0.000001
 * fm/c are representable internally.
 *
 * Usage:
 * ------
 * \code
 *   Clock labtime(0., 0.1);
 *   Clock endtime(10., 0.);
 *   while (labtime < endtime) {
 *     // do something
 *     // adapt the timestep size to external circumstances:
 *     if (system_is_very_dense()) {
 *       labtime.set_timestep_duration(labtime.timestep_duration() / 2.);
 *     }
 *     if (system_is_very_dilute()) {
 *       labtime.set_timestep_duration(labtime.timestep_duration() * 2.);
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
   * The value 0.000001 is very well suited because
   * \li It should be \f$10^{-n},\;n\in\mathbb{N}\f$. That's because we want to
   * use it to convert user input/output and that's in decimal representation.
   * \li The floating-point representation of the constant should have a small
   * error. 0.000001 has the smallest error (i.e. 0.022 ulp) in the range
   * \f$1\leq n \leq 10\f$. The small error helps to convert the internal
   * integer representation as precise as possible to floating-point.
   */
  static constexpr double resolution = 0.000001;

 public:
  /// The type used for counting ticks/time.
  using Representation = std::int64_t;

 public:
  /// default initializer: Timestep size is set to 0!
  Clock() = default;
  /**
   * Initialize with base time and time step size.
   *
   * \param[in] time base time
   * \param[in] dt step size
   */
  Clock(const double time, const double dt)
      : timestep_duration_(convert(dt)), reset_time_(convert(time)) {
    if (dt < 0.) {
      throw std::range_error("No negative time increment allowed");
    }
  }
  /// \return the current time.
  double current_time() const {
    return convert(reset_time_ + timestep_duration_ * counter_);
  }
  /**
   * \return the time in the next tick.
   *
   * This function is needed, because current_time() + timestep_duration()
   * is not the same as the next tick (numerically; this is due to
   * floating point arithmetic).
   */
  double next_time() const {
    if (counter_ * timestep_duration_ >=
        std::numeric_limits<Representation>::max() - timestep_duration_) {
      throw std::overflow_error("Too many timesteps, clock overflow imminent");
    }
    return convert(reset_time_ + timestep_duration_ * (counter_ + 1));
  }
  /// \return the time step size.
  double timestep_duration() const { return convert(timestep_duration_); }
  /**
   * Sets the time step size (and resets the counter).
   *
   * \param[in] dt new time step size
   */
  void set_timestep_duration(const double dt) {
    if (dt < 0.) {
      throw std::range_error("No negative time increment allowed");
    }
    reset_time_ += timestep_duration_ * counter_;
    counter_ = 0;
    timestep_duration_ = convert(dt);
  }
  /**
   * Checks if a multiple of a given interval is reached within the
   * next tick.
   *
   * \param[in] interval The interval \f$t_i\f$ at which, for instance,
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
  bool multiple_is_in_next_tick(const double interval) const {
    if (interval < 0.) {
      throw std::range_error("Negative interval makes no sense for clock");
    }
    const Representation int_interval = convert(interval);
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
  /**
   * Returns the next multiple of a given interval.
   *
   * \param[in] interval the given interval
   *
   * \return The smallest multiple of \p interval that is larger than
   * the current time.
   */
  double next_multiple(const double interval) const {
    const Representation int_interval = convert(interval);
    const auto current = reset_time_ + timestep_duration_ * counter_;
    if (unlikely(current < 0)) {
      return convert(current / int_interval * int_interval);
    }
    return convert((current / int_interval + 1) * int_interval);
  }
  /**
   * Set the time step such that it ends on the next multiple of the interval.
   *
   * \param interval The given interval
   */
  void end_tick_on_multiple(const double interval) {
    const Representation int_interval = convert(interval);
    const auto current = reset_time_ + timestep_duration_ * counter_;
    reset_time_ = current;
    counter_ = 0;
    if (unlikely(current < 0)) {
      timestep_duration_ = current / int_interval * int_interval - current;
    } else {
      timestep_duration_ =
          (current / int_interval + 1) * int_interval - current;
    }
  }
  /**
   * Resets the time to a pre-defined value \p reset_time.
   *
   * This is the only way of turning the clock back. It is needed so
   * that the time can be adjusted after initialization (different
   * initial conditions may require different starting times).
   *
   * \param[in] reset_time New time
   */
  void reset(const double reset_time) {
    if (reset_time < current_time()) {
      logger<LogArea::Clock>().debug("Resetting clock from", current_time(),
                                     " fm/c to ", reset_time, " fm/c");
    }
    reset_time_ = convert(reset_time);
    counter_ = 0;
  }
  /**
   * Advances the clock by one tick (\f$\Delta t\f$).
   *
   * This operator is used as `++clock`. The operator `clock++` is not
   * implemented deliberately, because that requires a copy of the clock
   * being created.
   */
  Clock& operator++() {
    // guard against overflow:
    if (counter_ >= std::numeric_limits<Representation>::max() - 1) {
      throw std::overflow_error("Too many timesteps, clock overflow imminent");
    }
    ++counter_;
    return *this;
  }
  /**
   * Advances the clock by an arbitrary timestep (multiple of 0.000001 fm/c).
   *
   * \tparam T type of the timestep
   * \param[in] big_timestep Timestep by which the clock is advanced.
   * \note It uses a template parameter only for disambiguation with the
   * overload below.
   */
  template <typename T>
  typename std::enable_if<std::is_floating_point<T>::value, Clock&>::type
  operator+=(T big_timestep) {
    if (big_timestep < 0.) {
      throw std::range_error("The clock cannot be turned back.");
    }
    reset_time_ += convert(big_timestep);
    return *this;
  }
  /**
   * advances the clock by an arbitrary number of ticks.
   *
   * \param[in] advance_several_timesteps Number of the timesteps added
   *                                      to the clock
   * \throw OverflowError if the number of the added timesteps exceeds
   *                      the maximum value.
   */
  Clock& operator+=(Representation advance_several_timesteps) {
    if (counter_ >= std::numeric_limits<Representation>::max() -
                        advance_several_timesteps) {
      throw std::overflow_error("Too many timesteps, clock overflow imminent");
    }
    counter_ += advance_several_timesteps;
    return *this;
  }
  /**
   * Compares the times between two clocks.
   *
   * \param[in] rhs The other clock.
   */
  bool operator<(const Clock& rhs) const {
    return (reset_time_ + timestep_duration_ * counter_) <
           (rhs.reset_time_ + rhs.timestep_duration_ * rhs.counter_);
  }
  /**
   * Compares the time of the clock against a fixed time.
   *
   * \param[in] time The other time.
   */
  bool operator<(double time) const { return current_time() < time; }
  /**
   * Compares the time of the clock against a fixed time.
   *
   * \param[in] time The other time.
   */
  bool operator>(double time) const { return current_time() > time; }

 private:
  /// A multiplier transfering the internal integer to the real time.
  static constexpr double to_double = resolution;
  /// A multiplier transfering the real time to the internal integer.
  static constexpr double from_double = 1. / resolution;

  /// Convert a double \p x into the internal int representation.
  static Representation convert(double x) {
    return std::round(x * from_double);
  }
  /// Convert an internal int value \p x into the double representation.
  static double convert(Representation x) { return x * to_double; }

  /**
   * Clock tick. This is purely internal and will be reset when the
   * timestep size is changed.
   */
  Representation counter_ = 0;
  /// The time step size \f$\Delta t\f$ in $10^{-3}$ fm.
  Representation timestep_duration_ = 0u;
  /// The time of last reset (when counter_ was set to 0).
  Representation reset_time_ = 0;
};

}  // namespace smash

#endif  // SRC_INCLUDE_CLOCK_H_
