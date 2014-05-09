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

TEST(compare) {
  Clock labtime(0.0f, 0.0f);
  Clock comtime(1.0f, 0.0f);
  VERIFY(labtime < comtime);
  VERIFY(labtime < 0.1f);
  VERIFY(comtime > 0.1f);
}

TEST(multiple) {
  Clock labtime(7.0f, 1.0f);
  float interval = 2.1f;
  // time is 7, next multiple is 8.4
  VERIFY(!labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime;  // 8; 8.4
  VERIFY( labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime;  // 9; 10.5
  VERIFY(!labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime;  // 10; 10.5
  VERIFY( labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime;  // 11; 12.6
  VERIFY(!labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime;  // 12; 12.6
  VERIFY( labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime;  // 13
  ++labtime;  // 14
  ++labtime;  // 15
  ++labtime;  // 16
  ++labtime;  // 17
  ++labtime;  // 18; 18.9
  VERIFY( labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime;  // 19; 21.0
  VERIFY(!labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime;  // 20; 21.0 (NO!)
  VERIFY(!labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime;  // 21; 21.0 (YES!)
  VERIFY( labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime;  // 22; 23.1 (NO!)
  VERIFY(!labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
}

TEST(multiple_small_interval) {
  Clock labtime(6.5f, 1.0f);
  float interval = 0.1f;
  for (int ticks = 0; ticks < 10; ++ticks) {
    VERIFY(labtime.multiple_is_in_next_tick(interval)) << ticks;
  }
}

TEST(assignment) {
  Clock labtime(4.2f, 0.3f);
  Clock resettime = labtime;
  ++labtime;
  FUZZY_COMPARE(labtime.current_time(), 4.5f);
  labtime = std::move(resettime);
  FUZZY_COMPARE(labtime.current_time(), 4.2f);
}
