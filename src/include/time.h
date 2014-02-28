/*
 *    Copyright (c) 2012-2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_TIME_H_
#define SRC_INCLUDE_TIME_H_

#include <time.h>

namespace Smash {

int clock_gettime(struct timespec *time);

}  // namespace Smash

#endif  // SRC_INCLUDE_TIME_H_
