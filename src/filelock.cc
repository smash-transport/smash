/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/filelock.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

namespace smash {

FileLock::FileLock(const bf::path& path) : path_(path), acquired_(false) {}

bool FileLock::acquire() {
  if (acquired_) {
    throw std::runtime_error("FileLock \"" + path_.native() +
                             "\" was already acquired.");
  }
  // The following POSIX syscall fails if the file already exists and creates it
  // otherwise. This is atomic, so there cannot be a race.
  const int fd =
      open(path_.native().c_str(), O_EXCL | O_CREAT, S_IRUSR | S_IWUSR);
  if (fd < 0) {
    return false;
  }
  acquired_ = true;
  if (close(fd) != 0) {
    throw std::runtime_error("Could not close file lock.");
  }
  return true;
}

FileLock::~FileLock() {
  if (acquired_) {
    bf::remove(path_);
  }
}

}  // namespace smash
