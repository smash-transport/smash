/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/fpenvironment.h"

#if defined __SSE__
#include <xmmintrin.h>
#endif

#include <csignal>

#include "smash/logging.h"

namespace smash {

#if !defined _GNU_SOURCE && defined __SSE__ && !defined __clang__
bool enable_float_traps(int femask) {
  static_assert(FE_INVALID == 0x01,
                "incorrect assumption that FE_INVALID == 0x01");
  static_assert(FE_DIVBYZERO == 0x04,
                "incorrect assumption that FE_DIVBYZERO == 0x04");
  static_assert(FE_OVERFLOW == 0x08,
                "incorrect assumption that FE_OVERFLOW == 0x08");
  static_assert(FE_UNDERFLOW == 0x10,
                "incorrect assumption that FE_UNDERFLOW == 0x10");
  static_assert(FE_INEXACT == 0x20,
                "incorrect assumption that FE_INEXACT == 0x20");

  // only accept supported values
  femask &= FE_ALL_EXCEPT;

#if defined __GNUC__ /*for inline asm*/
  // get the current FPU control word
  uint16_t fpucw;
  asm volatile("fstcw %0" : "=m"(fpucw));

  // clear the bits where the FPU should trap
  fpucw &= ~femask;

  // load the FPU control word back
  asm volatile("fldcw %0" ::"m"(fpucw));
#endif

  // the SSE CSR has the bit positions 7 bit positions further to the left
  femask <<= 7;

  /* Get the current CSR, mask of the relevant bits, and load it back into the
   * SSE CSR.*/
  _mm_setcsr(_mm_getcsr() & ~femask);

  // we did something, so return true
  return true;
}
#endif

void DisableFloatTraps::reenable_traps(int mask) {
  if (!enable_float_traps(mask)) {
    const auto &log = logger<LogArea::Fpe>();
    log.warn("Failed to setup traps on ", mask);
  }
}

void setup_default_float_traps() {
  {
    const auto &log = logger<LogArea::Fpe>();

    // pole error occurred in a floating-point operation:
    if (!enable_float_traps(FE_DIVBYZERO)) {
      log.warn("Failed to setup trap on pole error.");
    }

    // domain error occurred in an earlier floating-point operation:
    if (!enable_float_traps(FE_INVALID)) {
      log.warn("Failed to setup trap on domain error.");
    }

    /* The result of the earlier floating-point operation was too large to be
     * representable: */
    if (!enable_float_traps(FE_OVERFLOW)) {
      log.warn("Failed to setup trap on overflow.");
    }

    /* There's also FE_UNDERFLOW, where the result of the earlier
     * floating-point operation was subnormal with a loss of precision.
     * We do not consider this an error by default.
     * Furthermore, there's FE_INEXACT, but this traps if "rounding was
     * necessary to store the result of an earlier floating-point
     * operation". This is common and not really an error condition. */
  }

// Install the signal handler if we have the functionality.
#if (defined _POSIX_C_SOURCE && _POSIX_C_SOURCE >= 199309L) || \
    (defined _XOPEN_SOURCE && _XOPEN_SOURCE) ||                \
    (defined _POSIX_SOURCE && _POSIX_SOURCE)
/* The missing-fields-initializers warning says that not all fields of this
 * struct were explicitly initialized. That's exactly what is the intention
 * here, because then they are default-initialized, which is zero. */
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
  struct sigaction action = {};
  action.sa_flags = SA_SIGINFO;
  action.sa_sigaction = [](int signal, siginfo_t *info, void *) {
    const auto &log = logger<LogArea::Fpe>();
    if (signal == SIGFPE) {
      const char *msg = nullptr;
      switch (info->si_code) {
        case FPE_FLTDIV:
          msg = "Division by Zero (NaN)";
          break;
        case FPE_FLTUND:
          msg = "Underflow (result was subnormal with a loss of precision)";
          break;
        case FPE_FLTOVF:
          msg = "Overflow (result was too large to be representable)";
          break;
        case FPE_FLTINV:
          msg = "Invalid (domain error occurred)";
          break;
        case FPE_FLTRES:
          msg =
              "Inexact Result (rounding was necessary to store the result of "
              "an earlier floating-point operation)";
          break;
        default:
          msg = "unknown";
          break;
      }
      log.fatal("Floating point trap was raised: ", msg);
    } else {
      log.fatal("Unexpected Signal ", signal,
                " received in the FPE signal handler. Aborting.");
    }
    abort();
  };
  sigaction(SIGFPE, &action, nullptr);
#endif
}

}  // namespace smash
