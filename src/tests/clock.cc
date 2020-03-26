/*
 *
 *    Copyright (c) 2014-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <vir/test.h>  // This include has to be first

#include "../include/smash/clock.h"

using namespace smash;

TEST(size) {
  // if this fails, then either the internal structure in Clock is
  // changed (using other types or the addition of new variables) or the
  // alignment of the internal structure is somehow different. In both
  // cases, this test is meant to warn future developers to check if the
  // new behaviour is wanted or an unintendet side effect that should be
  // corrected.
  COMPARE(sizeof(Clock), 2 * sizeof(Clock::Representation));
  COMPARE(sizeof(UniformClock), 4 * sizeof(Clock::Representation));
}

TEST(set_clock) {
  UniformClock labtime(0.123, 0.234);
  COMPARE(labtime.current_time(), 0.123);
  FUZZY_COMPARE(labtime.timestep_duration(), 0.234);
}

TEST(run_clock) {
  UniformClock labtime(0.0, 0.1);
  COMPARE(labtime.current_time(), 0.0);
  ++labtime;
  FUZZY_COMPARE(labtime.current_time(), 0.1);
  labtime += 0.5;
  FUZZY_COMPARE(labtime.current_time(), 0.6);
  labtime += 2;
  FUZZY_COMPARE(labtime.current_time(), 0.8);
  UniformClock endtime(1.0, 0.0);
  while (labtime < endtime) {
    ++labtime;
  }
  FUZZY_COMPARE(labtime.current_time(), 1.0);
}

TEST(reset_timestep) {
  UniformClock labtime(0.0, 0.1);
  ++labtime;
  ++labtime;
  labtime.set_timestep_duration(0.2);
  ++labtime;
  FUZZY_COMPARE(labtime.current_time(), 0.4);
}

TEST(compare) {
  UniformClock labtime(0.0, 0.0);
  UniformClock comtime(1.0, 0.0);
  VERIFY(labtime < comtime);
  VERIFY(labtime < 0.1);
  VERIFY(comtime > 0.1);
}

TEST(assignment) {
  UniformClock labtime(4.2, 0.3);
  UniformClock resettime = labtime;
  ++labtime;
  COMPARE(labtime.current_time(), 4.5);
  labtime = std::move(resettime);
  COMPARE(labtime.current_time(), 4.2);
}

TEST_CATCH(init_negative_dt, std::range_error) {
  UniformClock labtime(4.4, -0.2);
}
TEST_CATCH(set_negative_dt, std::range_error) {
  UniformClock labtime(4.4, 0.2);
  labtime.set_timestep_duration(-0.3);
}
TEST_CATCH(big_timestep_negative, std::range_error) {
  UniformClock labtime(4.4, 0.2);
  labtime += -0.8;
}

TEST(no_overflow_single_increment) {
  UniformClock labtime(0.0, 1.0);
  labtime += (std::numeric_limits<Clock::Representation>::max() - 3);
  ++labtime;
  ++labtime;
}

TEST_CATCH(overflow_single_increment, std::overflow_error) {
  UniformClock labtime(0.0, 1.0);
  labtime += (std::numeric_limits<Clock::Representation>::max() - 3);
  ++labtime;
  ++labtime;
  ++labtime;
}

TEST(no_overflow_large_increment) {
  UniformClock labtime(0.0, 1.0);
  ++labtime;
  ++labtime;
  labtime += (std::numeric_limits<Clock::Representation>::max() - 3);
}

TEST_CATCH(overflow_large_increment, std::overflow_error) {
  UniformClock labtime(0.0, 1.0);
  ++labtime;
  ++labtime;
  ++labtime;
  labtime += (std::numeric_limits<Clock::Representation>::max() - 3);
}
