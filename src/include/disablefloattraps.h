/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_DISABLEFPE_H_
#define SRC_INCLUDE_DISABLEFPE_H_

#include <cfenv>
#include "logging.h"

namespace Smash {

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
   * \param reenable A bitwise or of the traps you want to keep enabled.
   */
  DisableFloatTraps(int reenable = 0) {
    std::feholdexcept(&environment_);
    if (reenable != 0) {
      if (-1 == feenableexcept(reenable)) {
        const auto &log = logger<LogArea::Fpe>();
        log.warn("Failed to setup traps on ", reenable);
      }
    }
  }

  /// When the guard goes out of scope the floating point environment is
  /// restored.
  ~DisableFloatTraps() { std::fesetenv(&environment_); }

 private:
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
 */
template <typename F>
void without_float_traps(F &&f) {
  DisableFloatTraps guard;
  f();
}

}  // namespace Smash

#endif  // SRC_INCLUDE_DISABLEFPE_H_
