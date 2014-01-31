/*
 *    Copyright (c) 2013
 *      SMASH Team
 * 
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <time.h>

#ifdef __MACH__
#include <sys/time.h>
#endif

#include "include/time.h"

#ifdef __MACH__
/* clock_gettime is not implemented use gettimeofday() */
int clock_gettime(struct timespec* time) {
  struct timeval now;
  int rv = gettimeofday(&now, NULL);
  /* return failure */
  if (rv)
    return rv;
  time->tv_sec  = now.tv_sec;
  time->tv_nsec = now.tv_usec * 1000;
  return 0;
}
#else
/* POSIX Linux and BSD clock_gettime() */
int clock_gettime(struct timespec* time) {
  return clock_gettime(CLOCK_REALTIME, time);
}
#endif
