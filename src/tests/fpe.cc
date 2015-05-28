/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "../include/disablefloattraps.h"

#include <csignal>

float blackhole = 0.f;
float divisor = 0.f;

static bool got_fpe = false;
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
  float a = 1.f;

  // default is no trap
  blackhole = a / divisor;
  COMPARE(got_fpe, false);

  // now it must trap
  feenableexcept(FE_DIVBYZERO);
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
