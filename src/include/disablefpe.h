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

class DisableFpe {
 public:
  DisableFpe(int reenable = 0) { std::feholdexcept(&environment_);
    if (reenable != 0) {
      if (-1 == feenableexcept(reenable)) {
        const auto &log = logger<LogArea::Fpe>();
        log.warn("Failed to setup traps on ", reenable);
      }
    }
  }

  ~DisableFpe() {
    std::fesetenv(&environment_);
  }

 private:
  fenv_t environment_;
};

template <typename F>
void float_environment_no_traps(F &&f) {
  DisableFpe guard;
  f();
}

}  // namespace Smash

#endif  // SRC_INCLUDE_DISABLEFPE_H_
