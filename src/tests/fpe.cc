/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "../include/fpenvironment.h"

#include <atomic>
#include <csignal>

float blackhole = 0.f;
float divisor = 0.f;

// this is atomic to enable the signal handler to change the value without the
// compiler assuming it cannot change between the two COMPAREs.
std::atomic<bool> got_fpe;

static void handle_fpe(int s) {
  if (s == SIGFPE) {
    got_fpe = true;
    std::feclearexcept(FE_ALL_EXCEPT);
    divisor = 1.f;
  }
}

TEST(setup_signal_handler) {
  std::signal(SIGFPE, &handle_fpe);
  COMPARE(got_fpe, false);
}

TEST(div_by_zero) {
  got_fpe = false;

  // default is no trap
  blackhole = 1.f / divisor;
  COMPARE(got_fpe, false);

  // now it must trap
  Smash::enable_float_traps(FE_DIVBYZERO);
  blackhole = 2.f / divisor;
  COMPARE(got_fpe, true);
  got_fpe = false;
  divisor = 0.f;

  // temporarily disable the trap
  Smash::without_float_traps([&] {
    blackhole = 3.f / divisor;
    COMPARE(got_fpe, false);
  });

  // and now trapping must work again
  blackhole = 4.f / divisor;
  COMPARE(got_fpe, true);
}
