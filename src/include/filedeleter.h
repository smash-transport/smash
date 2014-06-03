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

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <memory>
#include <stdexcept>

namespace Smash {
struct FileDeleter {
  constexpr FileDeleter() = default;
  void operator()(std::FILE *f) const {
    if (f == nullptr) {
      return;
    }
    if (0 != std::fclose(f)) {
      throw std::runtime_error(std::strerror(errno));
    }
  }
};

using FilePtr = std::unique_ptr<std::FILE, FileDeleter>;
}  // namespace Smash

#endif  // SRC_INCLUDE_FILEDELETER_H_
