/*
 *
 *    Copyright (c) 2014-2015,2017-2020,2022-2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/clock.h"

using namespace smash;

TEST(size) {
  // if this fails, then either the internal structure in Clock is
  // changed (using other types or the addition of new variables) or the
  // alignment of the internal structure is somehow different. In both
  // cases, this test is meant to warn future developers to check if the
  // new behaviour is wanted or an unintended side effect that should be
  // corrected.
  COMPARE(sizeof(Clock), 2 * sizeof(Clock::Representation));
  COMPARE(sizeof(UniformClock), 5 * sizeof(Clock::Representation));
}

TEST(set_clock) {
  UniformClock labtime(0.123, 0.234, 300.0);
  COMPARE(labtime.current_time(), 0.123);
  FUZZY_COMPARE(labtime.timestep_duration(), 0.234);
}

TEST(set_clock_negative_start_end_time) {
  UniformClock labtime(-0.123, 0.04, -0.03);
  COMPARE(labtime.current_time(), -0.123);
  FUZZY_COMPARE(labtime.timestep_duration(), 0.04);
}

TEST(run_clock) {
  UniformClock labtime(0.0, 0.1, 300.0);
  COMPARE(labtime.current_time(), 0.0);
  ++labtime;
  FUZZY_COMPARE(labtime.current_time(), 0.1);
  labtime += 0.5;
  FUZZY_COMPARE(labtime.current_time(), 0.6);
  labtime += 2;
  FUZZY_COMPARE(labtime.current_time(), 0.8);
  const double endtime = 1.0;
  while (labtime < endtime) {
    ++labtime;
  }
  FUZZY_COMPARE(labtime.current_time(), 1.0);
}

TEST(tick_clock_beyond_end_time) {
  const auto end_time = 10.0;
  UniformClock labtime(0.0, 10, end_time);
  VERIFY(labtime < end_time);
  ++labtime;
  VERIFY(!(labtime < end_time));
  VERIFY(!(labtime > end_time));
  ++labtime;
  VERIFY(labtime > end_time);
}

TEST(run_clock_across_end_time) {
  const auto end_time = 20.100001;
  UniformClock labtime(0.0, 0.1, end_time);
  auto counter = 0u;
  while (labtime < end_time) {
    ++labtime;
    ++counter;
  }
  COMPARE(counter, 202);
}

TEST(reset_timestep) {
  UniformClock labtime(0.0, 0.1, 0.6);
  ++labtime;
  ++labtime;
  labtime.set_timestep_duration(0.3);
  ++labtime;
  FUZZY_COMPARE(labtime.current_time(), 0.5);
  labtime.reset(-0.75, true);
  FUZZY_COMPARE(labtime.current_time(), -0.9);
}

TEST(compare) {
  UniformClock labtime(0.0, 0.1, 300.0);
  UniformClock comtime(1.0, 0.1, 300.0);
  VERIFY(labtime < comtime);
  VERIFY(labtime < 0.1);
  VERIFY(comtime > 0.1);
}

TEST(assignment) {
  UniformClock labtime(4.2, 0.3, 300.0);
  UniformClock resettime = labtime;
  ++labtime;
  COMPARE(labtime.current_time(), 4.5);
  labtime = std::move(resettime);
  COMPARE(labtime.current_time(), 4.2);
}

TEST_CATCH(create_clock_with_same_start_and_end_time, std::range_error) {
  UniformClock labtime(0.0, 10, 0.0);
}

TEST_CATCH(reset_time_bigger_than_end_time, std::range_error) {
  UniformClock labtime(4.4, 0.1, 1.0);
}

TEST_CATCH(init_zero_dt, std::range_error) {
  UniformClock labtime(4.4, 0.0, 300);
}

TEST_CATCH(init_negative_dt, std::range_error) {
  UniformClock labtime(4.4, -0.2, 300.0);
}

TEST_CATCH(set_negative_dt, std::range_error) {
  UniformClock labtime(4.4, 0.2, 300.0);
  labtime.set_timestep_duration(-0.3);
}

TEST_CATCH(big_timestep_negative, std::range_error) {
  UniformClock labtime(4.4, 0.2, 300.0);
  labtime += -0.8;
}

TEST(no_overflow_single_increment) {
  UniformClock labtime(0.0, 1.0, 300.0);
  labtime += (std::numeric_limits<Clock::Representation>::max() - 3);
  ++labtime;
  ++labtime;
}

TEST_CATCH(overflow_single_increment, std::overflow_error) {
  UniformClock labtime(0.0, 1.0, 300.0);
  labtime += (std::numeric_limits<Clock::Representation>::max() - 3);
  ++labtime;
  ++labtime;
  ++labtime;
}

TEST(no_overflow_large_increment) {
  UniformClock labtime(0.0, 1.0, 300.0);
  ++labtime;
  ++labtime;
  labtime += (std::numeric_limits<Clock::Representation>::max() - 3);
}

TEST_CATCH(overflow_large_increment, std::overflow_error) {
  UniformClock labtime(0.0, 1.0, 300.0);
  ++labtime;
  ++labtime;
  ++labtime;
  labtime += (std::numeric_limits<Clock::Representation>::max() - 3);
}

TEST_CATCH(custom_clock_tick_beyond_last, std::out_of_range) {
  CustomClock clock{{1, 2}};
  // Start time is 0.0 and not ticking the clock gives it
  VERIFY(clock.current_time() == 0);
  ++clock;
  VERIFY(clock.current_time() == 1);
  ++clock;
  VERIFY(clock.current_time() == 2);
  ++clock;
  clock.current_time();
}
TEST_CATCH(custom_clock_tick_till_end_ask_next_time, std::out_of_range) {
  CustomClock clock{{42}};
  ++clock;
  clock.next_time();
}

TEST_CATCH(custom_clock_remove_times_in_the_past, std::out_of_range) {
  CustomClock clock{{
      -3.14,
      -0.1,
      0.0,
      42,
  }};
  clock.remove_times_in_past(0.0);
  // Start time is 0.0 and not ticking the clock gives it
  VERIFY(clock.current_time() == 0.0);
  ++clock;
  VERIFY(clock.current_time() == 42);
  ++clock;
  clock.current_time();
}
