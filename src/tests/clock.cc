/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "tests/unittest.h"
#include "include/clock.h"

using namespace Smash;

TEST(set_clock) {
  Clock labtime(0.123f, 0.234f);
  COMPARE(labtime.current_time(), 0.123f);
  COMPARE(labtime.timestep_size(), 0.234f);
}

TEST(run_clock) {
  Clock labtime(0.0f, 0.1f);
  FUZZY_COMPARE(labtime.current_time(), 0.0f);
  ++labtime;
  FUZZY_COMPARE(labtime.current_time(), 0.1f);
  labtime += 0.5f;
  FUZZY_COMPARE(labtime.current_time(), 0.6f);
  labtime += 2;
  FUZZY_COMPARE(labtime.current_time(), 0.8f);
  Clock endtime(1.0f, 0.0f);
  while (labtime < endtime) {
    ++labtime;
  }
  FUZZY_COMPARE(labtime.current_time(), 1.0f);
}

TEST(reset_timestep) {
  Clock labtime(0.0f, 0.1f);
  ++labtime;
  ++labtime;
  labtime.set_timestep_size(0.2f);
  ++labtime;
  FUZZY_COMPARE(labtime.current_time(), 0.4f);
}
