/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

// FPE does not work with clang, so don't test it
#if !defined __clang__
#include "../include/smash/fpenvironment.h"

#include <csetjmp>
#include <csignal>
#include <stdexcept>

double blackhole = 0.;
double divisor = 0.;

std::jmp_buf jump_buffer;

struct fpe_exception : public std::exception {
  const char *what() const noexcept override {
    return "Floating Point Exception";
  }
};

static void handle_fpe(int s) {
  if (s == SIGFPE) {
    std::longjmp(jump_buffer, 1);
  }
}

static void do_division(double x) {
  if (setjmp(jump_buffer) == 0) {
    // normally goes here
    blackhole = x / divisor;
  } else {
    // longjmp goes here
    throw fpe_exception();
  }
}

// this test cannot be merged with the fpe.cc test because a single process can
// only handle a single SIGFPE. Actually program resumption after the SIGFPE
// handler is undefined behavior. But we're still in compiler extension land, so
// the test works.
TEST(without_float_traps) {
  std::signal(SIGFPE, &handle_fpe);

  smash::enable_float_traps(FE_DIVBYZERO);     // now div by zero must trap
  smash::without_float_traps([&] {             // temporarily disable the trap
    VERIFY(!std::fetestexcept(FE_DIVBYZERO));  // flag not set yet
    do_division(2.);  // this sets the flag, but doesn't trap
    VERIFY(std::fetestexcept(FE_DIVBYZERO));  // flag must be set now
  });
  VERIFY(!std::fetestexcept(
      FE_DIVBYZERO));  // after the lambda the flag must be clear again

  bool got_exception = false;
  try {
    do_division(3.);  // after the lambda this must trap again
    FAIL() << "may never be reached";
  } catch (fpe_exception &) {  // the signal handler produces an exception via
                               // longjmp
    got_exception = true;
  }
  VERIFY(!std::fetestexcept(
      FE_DIVBYZERO));  // flag must not be set because it trapped
  VERIFY(got_exception);
}
#else
TEST(unsupported) {}
#endif
