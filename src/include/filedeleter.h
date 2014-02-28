/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_FILEDELETER_H_
#define SRC_INCLUDE_FILEDELETER_H_

#include <cstdio>
#include <cstring>
#include <memory>
#include <stdexcept>

namespace std {
template <>
struct default_delete<std::FILE> {
  constexpr default_delete() = default;
  void operator()(std::FILE *f) const {
    if (f == nullptr) {
      return;
    }
    if (0 != std::fclose(f)) {
      throw runtime_error(strerror(errno));
    }
  }
};
}  // namespace std

#endif  // SRC_INCLUDE_FILEDELETER_H_
