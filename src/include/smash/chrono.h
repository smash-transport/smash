/*
 *    Copyright (c) 2013-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_CHRONO_H_
#define SRC_INCLUDE_CHRONO_H_

#include <chrono>

/**
 * \file
 *
 * Collection of useful type aliases to measure and output the (real)
 * runtime. Not connected to the simulation time.
 */

namespace smash {

/// Type (alias) that is used to store the current time.
using SystemTimePoint = std::chrono::time_point<std::chrono::system_clock>;

/// Type (alias) used to obtain the current time via SystemClock:Now().
using SystemClock = std::chrono::system_clock;

/// The time duration type (alias) used for measuring run times.
using SystemTimeSpan = SystemClock::duration;

}  // namespace smash

#endif  // SRC_INCLUDE_CHRONO_H_
