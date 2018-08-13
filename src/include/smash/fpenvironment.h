/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_FPENVIRONMENT_H_
#define SRC_INCLUDE_FPENVIRONMENT_H_

#include <cfenv>

namespace smash {

#if defined _GNU_SOURCE && !defined __clang__
/**
 * Standard C/C++ don't have a function to modify the trapping behavior. You
 * can only save and restore the setup. With glibc you can change it via
 * feenableexcept and fedisableexcept. Without glibc inline asm and SSE
 * intrinsics can do it (for x86).
 *
 * It is important that the used compiler supports floating-point traps, other-
 * wise it might break them by reordering floating-point operations when
 * optimizing. Currently floating-point traps are disabled for Clang, because
 * it does not support them.
 *
 * \param mask A bitwise of the traps you want to keep enabled.
 * \return Whether the trap was successfully set.
 *
 * glibc specific implementation
 */
inline bool enable_float_traps(int mask) { return -1 != feenableexcept(mask); }
#elif defined __SSE__ && !defined __clang__
/// Directly program the trap on the SSE unit
bool enable_float_traps(int mask);
#else
/// Fallback that fails to set the trap
inline bool enable_float_traps(int) { return false; }
#endif

/**
 * Setup the floating-point traps used throughout SMASH.
 *
 * If possible, this function additionally installs a signal handler that prints
 * what kind of condition triggered the trap. This requires POSIX.1-2001 to
 * work.
 */
void setup_default_float_traps();

/**
 * Guard type that safely disables floating point traps for the scope in which
 * it is placed.
 *
 * Example:
 * \code
 * {
 *   // some code where FPEs will trap
 *   DisableFloatTraps guard;
 *   // all code up to the closing brace will not trap anymore.
 * }
 * \endcode
 *
 * You can also keep some traps enabled. E.g. FE_DIVBYZERO is a candidate you
 * might want to keep enabled, whereas FE_UNDERFLOW and FE_OVERFLOW are the ones
 * that you really need to get rid of:
 * \code
 * DisableFloatTraps guard(FE_DIVBYZERO | FE_INVALID);
 * \endcode
 *
 * \note There is no guarantee about the complexity of modifying the floating
 * point environment. This could be very expensive. Therefore it is always
 * preferable to fix your code to not require underflow or overflow.
 *
 * \see http://en.cppreference.com/w/cpp/numeric/fenv/FE_exceptions for a
 * complete list of flags.
 */
class DisableFloatTraps {
 public:
  /**
   * Constructs the guard object.
   *
   * \param mask A bitwise of the traps you want to keep enabled.
   * \return The constructed object.
   */
  explicit DisableFloatTraps(int mask = 0) {
    std::feholdexcept(&environment_);
    if (mask != 0) {
      reenable_traps(mask);
    }
  }

  /**
   * When the guard goes out of scope the floating point environment is
   * restored.
   */
  ~DisableFloatTraps() { std::fesetenv(&environment_); }

 private:
  /**
   * Reenables the given traps.
   * \param mask A bitwise of the traps you want to keep enabled.
   */
  void reenable_traps(int mask);

  /// The stored environment that the destructor will restore.
  std::fenv_t environment_;
};

/**
 * Convenience function to create a scope where all floating point traps are
 * disabled.
 *
 * Example:
 * \code
 * // some code where FPEs will trap
 * without_float_traps([&] {
 *   // all code up to the closing brace will not trap anymore.
 * });
 * // more code where FPEs will trap
 * \endcode
 *
 * \param f A functor (e.g. lambda) that is executed in the cleared floating
 * point environment.
 * \tparam F type of the functor
 */
template <typename F>
void without_float_traps(F &&f) {
  DisableFloatTraps guard;
  f();
}

}  // namespace smash

#endif  // SRC_INCLUDE_FPENVIRONMENT_H_
