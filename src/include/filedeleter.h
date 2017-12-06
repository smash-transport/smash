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

namespace smash {

/**
 * FileDeleter is the deleter class for std::unique_ptr of std::FILE.
 *
 * std::unique_ptr takes a second template argument which determines what
 * happens when the resource it holds needs to be freed. The default
 * implementation calls `delete`. For std::FILE the resource needs to be freed
 * with a call to std::fclose instead. Therefore FilePtr requires a custom
 * deleter class to correctly free the resource.
 */
struct FileDeleter {
  /// The class has no members, so this is a noop.
  constexpr FileDeleter() = default;

  /// frees the std::FILE resource if it is non-zero.
  void operator()(std::FILE *f) const {
    if (f == nullptr) {
      return;
    }
    if (0 != std::fclose(f)) {
      throw std::runtime_error(std::strerror(errno));
    }
  }
};

/**
 * A RAII type to replace `std::FILE *`.
 *
 * This is an alias type for std::unique_ptr to automatically free the std::FILE
 * resource after the last reference goes out of scope. It is important to use a
 * custom deleter type, and therefore SMASH code should never use
 * std::unique_ptr directly with std::FILE.
 */
using FilePtr = std::unique_ptr<std::FILE, FileDeleter>;
}  // namespace smash

#endif  // SRC_INCLUDE_FILEDELETER_H_
