/*
 *
 *    Copyright (c) 2014-2020,2022-2025
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
#include "numeric_cast.h"

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
 * \attention
 * There are two very different types of clocks implemented.
 * -# The UniformClock is meant to track time in a simulation and it has the
 *    quite peculiar feature that it is aware of the end time of the simulation.
 *    Therefore, although it can be ticked beyond such a time, its methods about
 *    the current time, the next time and the time-step duration will adjust
 *    their return value when the clock gets to or past the simulation end time.
 * -# The CustomClock is meant to track a list of time points and it is best
 *    suited for output. Ticking this clock means to consider the following
 *    point in time. Ticking beyond the last time point is considered an error.
 *
 * <h3> Potential usage for adapting time steps: </h3>
 *
 * \code
 *   const double end_time = 10.;
 *   UniformClock lab_time(0., 0.1, end_time);
 *   while (lab_time < end_time) {
 *     // do something
 *     // adapt the timestep size to external circumstances:
 *     if (system_is_very_dense()) {
 *       lab_time.set_timestep_duration(lab_time.timestep_duration() / 2.);
 *     }
 *     if (system_is_very_dilute()) {
 *       lab_time.set_timestep_duration(lab_time.timestep_duration() * 2.);
 *     }
 *     // let the clock tick
 *     ++lab_time;
 *   }
 * \endcode
 *
 * <h3> Possible actions for Clock </h3>
 *
 * \li look at it and find out the current time
 * \see current_time()
 * \see next_time()
 * \li advance the clock (by one tick, by several ticks, or by a given
 * time)
 * \see operator++()
 * \see UniformClock::operator+=(T)
 * \see operator+=(Representation)
 * \li set / retrieve the timestep (length of one tick)
 * \see UniformClock::set_timestep_duration(double)
 * \see timestep_duration()
 * \li compare time against different clock or fixed value
 * \see operator<(const Clock&) const
 * \see operator<(double) const
 * \see operator>(double) const
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
   * Reset the clock to the starting time of the simulation
   *
   * \param[in] start_time starting time of the simulation
   * \param[in] is_output_clock whether this is an output clock rather than a
   *                            lab clock
   */
  virtual void reset(double start_time, bool is_output_clock) = 0;
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
   * \param[in] advance_several_timesteps Number of the time steps added
   *                                      to the clock
   * \throw OverflowError if the number of the added time steps exceeds
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
   * Compares the internal times of two clocks.
   *
   * \param[in] rhs The other clock.
   */
  bool operator<(const Clock& rhs) const {
    return present_internal_time() < rhs.present_internal_time();
  }

  /**
   * Compares the internal time of the clock against a fixed time.
   *
   * \param[in] time The other time.
   */
  bool operator<(double time) const { return present_internal_time() < time; }

  /**
   * Compares the internal time of the clock against a fixed time.
   *
   * \param[in] time The other time.
   */
  bool operator>(double time) const { return present_internal_time() > time; }

  virtual ~Clock() = default;

 protected:
  /**
   * This function \b always returns the clock time, even if children might
   * attribute a different behaviour to \c current_time method (as UniformClock
   * does).
   *
   * \warning It is important to have this method that is in turn used in the
   * comparison operators, so that clock comparisons are independent from
   * other possible existing mechanism (like that of the UniformClock).
   *
   * \return The present internal clock time.
   */
  virtual double present_internal_time() const = 0;
  /**
   * Internally used to count the number of time steps.
   */
  Representation counter_ = 0;
};

/** Clock with uniformly spaced time steps
 *
 * <h3> Internal clock mechanisms </h3>
 *
 * This clock stores a time step size \f$\Delta t\f$, a base time \f$t_0\f$ as
 * well as an end time \f$t_{end}\f$ and a counter \f$n\f$. The current time is
 * calculated from \f$t = t_0 + n \cdot \Delta t\f$. When \f$\Delta t\f$ is
 * changed, \f$t_0\f$ is reset to the present time and \f$n\f$ is set to 0.
 * As soon as \f$t\geq t_{end}\f$, the clock will always return \f$t_{end}\f$ as
 * current and next time. The time step size can be retrieved and the returned
 * value is
 * \f[
 * \begin{aligned}
 *   \Delta t       \qquad&\mbox{if}& &t\leq t_{end}-\Delta t \\
 *   t_{end}-t      \qquad&\mbox{if}& t_{end}-\Delta t<{}&t<t_{end} \\
 *   0.0            \qquad&\mbox{if}& &t\geq t_{end} \\
 * \end{aligned}
 * \f]
 * In the last case, i.e. if the time step size is required when the clock
 * ticked beyond the simulation end, a warning is given to the user.
 */
class UniformClock : public Clock {
  /**
   * Defines the resolution of the clock (namely the smallest representable time
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
  UniformClock(double time, double dt, double time_end)
      : timestep_duration_(convert(dt)),
        reset_time_(convert(time)),
        time_end_(convert(time_end)) {
    if (dt <= 0.) {
      throw std::range_error("Time increment must be positive and non-zero");
    }
    if (reset_time_ >= time_end_) {
      throw std::range_error(
          "The initial time of UniformClock must be smaller than the end time. "
          "(Attempt to set initial time to " +
          std::to_string(time) + " and end time to " +
          std::to_string(time_end) + " not possible)");
    }
  }
  /**
   * \return the current time or the end time if the clock ticked beyond it.
   */
  double current_time() const override {
    auto present_time = present_internal_time();
    // Do comparison in internal representation unit and return converted values
    if (convert(present_time) > time_end_) {
      return convert(time_end_);
    } else {
      return present_time;
    }
  }
  /**
   * \return the time in the next tick or the end time if the clock ticked
   *         beyond it.
   *
   * \note This function is needed, because current_time() + timestep_duration()
   *       is not the same as the next tick (numerically; this is due to
   *       floating point arithmetic).
   */
  double next_time() const override {
    if (counter_ * timestep_duration_ >=
        std::numeric_limits<Representation>::max() - timestep_duration_) {
      throw std::overflow_error("Too many timesteps, clock overflow imminent");
    }
    auto next_point_in_time = reset_time_ + timestep_duration_ * (counter_ + 1);
    if (next_point_in_time > time_end_) {
      return convert(time_end_);
    } else {
      return convert(next_point_in_time);
    }
  }

  /**
   * \return the time step size from the current time. If a full tick would
   * result in a time larger then the end time, a smaller size is returned. If
   * the clock is already beyond the end time, 0.0 is returned and a warning is
   * printed. \see UniformClock description.
   */
  double timestep_duration() const override {
    auto present_time = convert(present_internal_time());
    if (present_time > time_end_) {
      logg[LClock].warn() << "UniformClock asked for timestep duration beyond "
                             "end of simulation, returning 0.";
      return 0.0;
    } else if (present_time + timestep_duration_ > time_end_) {
      return convert(time_end_ - present_time);
    } else {
      return convert(timestep_duration_);
    }
  }
  /**
   * Sets the time step size (and resets the counter).
   *
   * \param[in] dt new time step size
   */
  void set_timestep_duration(double dt) {
    if (dt <= 0.) {
      throw std::range_error("Time increment must be positive and non-zero!");
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
  void reset(double start_time, bool is_output_clock) override {
    double reset_time;
    if (is_output_clock) {
      auto delta_t = convert(timestep_duration_);
      reset_time = std::floor(start_time / delta_t) * delta_t;
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
   * \param[in] big_timestep Time step by which the clock is advanced.
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

 protected:
  /**
   * Access the internal time of the clock, independently from the end time.
   *
   * \return the internal clock time.
   */
  double present_internal_time() const override {
    return convert(reset_time_ + timestep_duration_ * counter_);
  }

 private:
  /// A multiplier transferring the internal integer to the real time.
  static constexpr double to_double = resolution;
  /// A multiplier transferring the real time to the internal integer.
  static constexpr double from_double = 1. / resolution;

  /// Convert a double \p x into the internal int representation.
  static Representation convert(double x) {
    return numeric_cast<Representation>(std::round(x * from_double));
  }
  /// Convert an internal int value \p x into the double representation.
  static double convert(Representation x) { return x * to_double; }

  /// The time step size \f$\Delta t\f$ in \f$10^{-6}\,\mathrm{fm}\f$.
  Representation timestep_duration_ = 0;
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
   * \return The start time if the clock has never been ticked or the current
   *         time otherwise.
   * \throw std::out_of_range if the clock has ticked beyond the last time.
   * \throw std::runtime_error if the clock has an internal broken state.
   */
  double current_time() const override {
    if (counter_ == -1) {
      // current time before the first output should be the starting time
      return start_time_;
    } else if (counter_ < -1) {
      throw std::runtime_error(
          "Trying to access time of clock in invalid state.");
    } else {
      return custom_times_.at(counter_);
    }
  }

  /**
   * \return The next custom time.
   * \throw std::out_of_range if the clock has ticked beyond last time.
   */
  double next_time() const override { return custom_times_.at(counter_ + 1); }

  /// \copydoc Clock::timestep_duration
  double timestep_duration() const override {
    return next_time() - current_time();
  }

  /**
   * Reset the clock to the starting time of the simulation.
   *
   * \param[in] start_time starting time of the simulation
   *
   * \note The second \c bool parameter is irrelevant and unused here.
   */
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
    custom_times_.erase(
        std::remove_if(
            custom_times_.begin(), custom_times_.end(),
            [start_time](double t) {
              if (t < start_time) {
                logg[LClock].warn("Removing custom output time ", t,
                                  " fm since it is earlier than the "
                                  "starting time of the simulation");
                return true;
              } else if (t == start_time) {
                logg[LClock].debug(
                    "The start time ", t,
                    " fm has to be removed from the 'custom_times_' vector "
                    "since it will be otherwise considered twice in the actual "
                    "output time steps");
                return true;
              } else {
                return false;
              }
            }),
        custom_times_.end());
  }

 protected:
  /**
   * For the CustomClock, the internal time is basically by design the same as
   * what the current_time() method returns.
   *
   * \return the same as \c current_time does.
   */
  double present_internal_time() const override { return current_time(); }

 private:
  /// Vector of times where output is generated
  std::vector<double> custom_times_;
  /// Starting time of the simulation
  double start_time_ = 0.;
};
}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_CLOCK_H_
