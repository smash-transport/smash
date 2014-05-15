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

TEST(size) {
  // if this fails, then either the internal structure in Clock is
  // changed (using other types or the addition of new variables) or the
  // alignment of the internal structure is somehow different. In both
  // cases, this test is meant to warn future developers to check if the
  // new behaviour is wanted or an unintendet side effect that should be
  // corrected.
  size_t floatsize = sizeof(float);
  size_t int32size = sizeof(uint32_t);
  COMPARE(sizeof(Clock), 2 * floatsize + int32size);
}

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
  labtime += 2u;
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
  for (int ticks = 0; ticks < 20; ++ticks) {
    VERIFY(labtime.multiple_is_in_next_tick(interval)) << ticks;
  }
}

TEST(multiple_small_increment) {
  Clock labtime(0.0f, 0.002f);
  for (int i = 0; i < 30; ++i) {
    while (++labtime < i + 0.997f) {
      VERIFY(!labtime.multiple_is_in_next_tick(1.0f)) << labtime.current_time();
    }
    // check that one of this or next tick (but one only!) fits.
    VERIFY((labtime.multiple_is_in_next_tick(1.0f)
       ^ (++labtime).multiple_is_in_next_tick(1.0f))) << labtime.current_time();
  }
}

TEST(multiple_negative_times) {
  Clock labtime(-12.0f, 1.0f);
  float interval = 2.2f;
  VERIFY(!labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime; // time is -11, next multiple is -11.0: YES
  VERIFY( labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime; // -10, -8.8: NO
  VERIFY(!labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime; // -9; -8.8: YES
  VERIFY( labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime; // -8; -6.6: NO
  VERIFY(!labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime; // -7; -6.6: YES
  VERIFY( labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime; // -6; -4.4: NO
  VERIFY(!labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime; // -5; -4.4: YES
  VERIFY( labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime; // -4; -2.2: NO
  VERIFY(!labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime; // -3; -2.2: YES
  VERIFY( labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime; // -2;  0.0: NO
  VERIFY(!labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime; // -1;  0.0: NO
  VERIFY(!labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime; // 0;  0.0: YES
  VERIFY( labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime; // 1;  2.2: NO
  VERIFY(!labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
  ++labtime; // 2;  2.2: YES
  VERIFY( labtime.multiple_is_in_next_tick(interval)) << labtime.current_time();
}

TEST(assignment) {
  Clock labtime(4.2f, 0.3f);
  Clock resettime = labtime;
  ++labtime;
  FUZZY_COMPARE(labtime.current_time(), 4.5f);
  labtime = std::move(resettime);
  FUZZY_COMPARE(labtime.current_time(), 4.2f);
}

TEST_CATCH(init_negative_dt, std::range_error) {
  Clock labtime(4.4f, -0.2f);
}
TEST_CATCH(set_negative_dt, std::range_error) {
  Clock labtime(4.4f, 0.2f);
  labtime.set_timestep_size(-0.3f);
}
TEST_CATCH(negative_interval, std::range_error) {
  Clock labtime(4.4f, 0.2f);
  labtime.multiple_is_in_next_tick(-0.8f);
}
TEST_CATCH(big_timestep_negative, std::range_error) {
  Clock labtime(4.4f, 0.2f);
  labtime += -0.8f;
}

TEST_CATCH(overflow_single_increment, std::overflow_error) {
  Clock labtime(0.0f, 1.0f);
  labtime += (std::numeric_limits<uint32_t>::max() - 3);
  ++labtime;
  ++labtime;
  ++labtime;
  ++labtime;
  ++labtime;
  ++labtime;
}

TEST_CATCH(overflow_large_increment, std::overflow_error) {
  Clock labtime(0.0f, 1.0f);
  ++labtime;
  ++labtime;
  ++labtime;
  ++labtime;
  labtime += (std::numeric_limits<uint32_t>::max() - 3);
}
