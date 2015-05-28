/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_FPENVIRONMENT_H_
#define SRC_INCLUDE_FPENVIRONMENT_H_

#include <cfenv>

namespace Smash {

/**
 * Standard C/C++ don't have a function to modify the trapping behavior. You
 * can only save and restore the setup. With glibc you can change it via
 * feenableexcept and fedisableexcept. Without glibc inline asm and SSE
 * intrinsics can do it (for x86).
 */
#if defined _GNU_SOURCE
// glibc specific implementation
inline bool enable_float_traps(int mask) {
  return -1 != feenableexcept(mask);
}
#elif defined __SSE__
// directly program the trap on the SSE unit
bool enable_float_traps(int femask);
#else
// fallback that fails to set the trap
inline bool enable_float_traps(int) { return false; }
#endif

}  // namespace Smash

#endif  // SRC_INCLUDE_FPENVIRONMENT_H_
