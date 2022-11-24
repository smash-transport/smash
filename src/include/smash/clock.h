/*
 *
 *    Copyright (c) 2014-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_CLOCK_H_
#define SRC_INCLUDE_SMASH_CLOCK_H_

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

#include "logging.h"

namespace smash {

static constexpr int LClock = LogArea::Clock::id;

/**
 * Clock tracks the time in the simulation.
 *
 * The basic unit is 1 fm in natural units which correspond to
 * \f$\frac{10^{-15}}{299\,792\,458}\,\mathrm{s} \approx
 * 0.33\cdot10^{-23}\,\mathrm{s}\f$ in the international system of units. The
 * resolution of the clock is \f$0.000001\,\mathrm{fm}=10^{-6}\,\mathrm{fm}\f$,
 * i.e. only multiples of \f$0.000001\,\mathrm{fm}\f$ are internally
 * representable.
 *
 * Potential usage for adapting time steps:
 * ------
 * \code
 *   UniformClock labtime(0., 0.1, end_time_);
 *   UniformClock endtime(10., 0., end_time_);
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
 **/
class Clock {
 public:
  /// The type used for counting ticks/time.
  using Representation = std::int64_t;
  /// \return the duration of the current time step
  virtual double timestep_duration() const = 0;
  /// \return the current time
  virtual double current_time() const = 0;
  /// \return the time of the next time step
  virtual double next_time() const = 0;
  /**
   * reset the clock to the starting time of the simulation
   *
   * \param[in] start_time starting time of the imulation
   * \param[in] is_output_clock whether this is an output clock rather than a
   *                            lab clock
   */
  virtual void reset(double start_time, const bool is_output_clock) = 0;
  /**
   * Remove output times before the starting time of the simulation if this
   * is a custom clock.
   *
   * \param[in] start_time starting time of the simulation
   */

  virtual void remove_times_in_past(double start_time) = 0;
  /**
   * Advances the clock by one tick.
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
   * Advances the clock by an arbitrary number of ticks.
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
    return current_time() < rhs.current_time();
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

  virtual ~Clock() = default;

 protected:
  /**
   * Internally used to count the number of time steps.
   */
  Representation counter_ = 0;
};

/** Clock with uniformly spaced time steps
 *
 * Internals
 * ---------
 *
 * Clock stores a time step size \f$\Delta t\f$ and a base time
 * \f$t_0\f$ as well as a counter \f$n\f$. The current time is
 * calculated from \f$t = t_0 + n \cdot \Delta t\f$. When \f$\Delta t\f$
 * is changed, \f$t_0\f$ is reset.
 *
 */
class UniformClock : public Clock {
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
  /// default initializer: Timestep size is set to 0!
  UniformClock() = default;
  /**
   * Initialize with base time and time step size.
   *
   * \param[in] time base time
   * \param[in] dt step size
   * \param[in] time_end end time of particle propagation
   */
  UniformClock(const double time, const double dt, const double time_end)
      : timestep_duration_(convert(dt)),
        reset_time_(convert(time)),
        time_end_(convert(time_end)) {
    if (dt < 0.) {
      throw std::range_error("No negative time increment allowed");
    }
  }
  /// \return the current time.
  double current_time() const override {
    return convert(reset_time_ + timestep_duration_ * counter_);
  }
  /**
   * \return the time in the next tick.
   *
   * This function is needed, because current_time() + timestep_duration()
   * is not the same as the next tick (numerically; this is due to
   * floating point arithmetic).
   */
  double next_time() const override {
    if (counter_ * timestep_duration_ >=
        std::numeric_limits<Representation>::max() - timestep_duration_) {
      throw std::overflow_error("Too many timesteps, clock overflow imminent");
    }
    if ((reset_time_ + timestep_duration_ * (counter_ + 1)) > time_end_) {
      return convert(time_end_);
    } else {
      return convert(reset_time_ + timestep_duration_ * (counter_ + 1));
    }
  }
  /// \return the time step size.
  double timestep_duration() const override {
    if ((reset_time_ + timestep_duration_ * (counter_ + 1)) > time_end_) {
      Representation last_timestep =
          time_end_ - (reset_time_ + timestep_duration_ * counter_);
      return convert(last_timestep);
    } else {
      return convert(timestep_duration_);
    }
  }
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
   * Resets the time to the starting time of an event.
   *
   * \param[in] start_time Starting time of the simulation
   * \param[in] is_output_clock whether this is an output clock or a lab clock
   */
  void reset(const double start_time, const bool is_output_clock) override {
    double reset_time;
    if (is_output_clock) {
      reset_time =
          std::floor(start_time / timestep_duration()) * timestep_duration();
    } else {
      reset_time = start_time;
    }
    if (reset_time < current_time()) {
      logg[LClock].debug("Resetting clock from", current_time(), " fm to ",
                         reset_time, " fm");
    }
    reset_time_ = convert(reset_time);
    counter_ = 0;
  }

  void remove_times_in_past(double) override{};

  /**
   * Advances the clock by an arbitrary timestep (multiple of 0.000001 fm).
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

  /// The time step size \f$\Delta t\f$ in \f$10^{-6}\,\mathrm{fm}\f$.
  Representation timestep_duration_ = 0u;
  /// The time of last reset (when counter_ was set to 0).
  Representation reset_time_ = 0;
  /// The end time of the particle propagation
  Representation time_end_ = 0;
};

/// Clock with explicitly defined time steps
class CustomClock : public Clock {
 public:
  /**
   * Initialises a custom clock with explicitly given output times
   *
   * \param[in] times vector of desired output times
   */
  explicit CustomClock(std::vector<double> times) : custom_times_(times) {
    std::sort(custom_times_.begin(), custom_times_.end());
    counter_ = -1;
  }
  /**
   * \copydoc Clock::current_time
   * \throw runtime_error if the clock has never been advanced
   */
  double current_time() const override {
    if (counter_ == -1) {
      // current time before the first output should be the starting time
      return start_time_;
    } else if (counter_ < -1) {
      throw std::runtime_error("Trying to access undefined zeroth output time");
    } else {
      return custom_times_[counter_];
    }
  }
  /// \copydoc Clock::next_time
  double next_time() const override { return custom_times_[counter_ + 1]; }
  double timestep_duration() const override {
    return next_time() - current_time();
  }
  void reset(double start_time, bool) override {
    counter_ = -1;
    start_time_ = start_time;
  }

  /**
   * Remove all custom times before start_time.
   *
   * \param[in] start_time starting time of the simulation
   */
  void remove_times_in_past(double start_time) override {
    std::remove_if(custom_times_.begin(), custom_times_.end(),
                   [start_time](double t) {
                     if (t <= start_time) {
                       logg[LClock].warn("Removing custom output time ", t,
                                         " fm since it is earlier than the "
                                         "starting time of the simulation");
                       return true;
                     } else {
                       return false;
                     }
                   });
  }

 private:
  /// Vector of times where output is generated
  std::vector<double> custom_times_;
  /// Starting time of the simulation
  double start_time_ = 0.;
};
}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_CLOCK_H_
