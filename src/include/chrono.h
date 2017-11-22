/*
 *    Copyright (c) 2013-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_CHRONO_H_
#define SRC_INCLUDE_CHRONO_H_

#include <chrono>

namespace smash {
using SystemTimePoint = std::chrono::time_point<std::chrono::system_clock>;
using SystemClock = std::chrono::system_clock;
/**
 * The time duration type used for measuring run times.
 */
using SystemTimeSpan = SystemClock::duration;
}

#endif  // SRC_INCLUDE_CHRONO_H_
