/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/fpenvironment.h"
#include "include/logging.h"

#if defined __SSE__
#include <xmmintrin.h>
#endif

namespace Smash {

#if !defined _GNU_SOURCE && defined __SSE__
bool enable_float_traps(int femask) {
  static_assert(FE_INVALID == 0x01, "incorrect assumption that FE_INVALID == 0x01");
  static_assert(FE_DIVBYZERO == 0x04, "incorrect assumption that FE_DIVBYZERO == 0x04");
  static_assert(FE_OVERFLOW == 0x08, "incorrect assumption that FE_OVERFLOW == 0x08");
  static_assert(FE_UNDERFLOW == 0x10, "incorrect assumption that FE_UNDERFLOW == 0x10");
  static_assert(FE_INEXACT == 0x20, "incorrect assumption that FE_INEXACT == 0x20");

  // only accept supported values
  femask &= FE_ALL_EXCEPT;

#if defined __GNUC__ /*for inline asm*/
  // get the current FPU control word
  unsigned short fpucw;
  asm volatile("fstcw %0" : "=m"(fpucw));

  // clear the bits where the FPU should trap
  fpucw &= ~femask;

  // load the FPU control word back
  asm volatile("fldcw %0" :: "m"(fpucw));
#endif

  // the SSE CSR has the bit positions 7 bit positions further to the left
  femask <<= 7;

  // get the current CSR, mask of the relevant bits, and load it back into the
  // SSE CSR
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

}  // namespace Smash
