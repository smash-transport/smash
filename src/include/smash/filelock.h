/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_FILELOCK_H_
#define SRC_INCLUDE_FILELOCK_H_

#include <boost/filesystem.hpp>

#include "forwarddeclarations.h"

namespace smash {

/** Guard to create a file lock.
 *
 * This can be used to make sure concurrent processes do not interfere when for
 * example writing files.
 *
 * The lock file will be deleted when the FileLock is destroyed.
 *
 * Internally this uses some POSIX system calls to atomically check for the
 * existence of a file and create it if it does not exists. This works even
 * over NFS (for at least NFSv3 and Linux 2.6).
 */
class FileLock {
 public:
  /** Construct a file lock guard with a lock file at the given path.
   *
   * This will not create a file, use acquire() for that.
   *
   * \param[in] path Path to file.
   * \return Constructed object.
   */
  explicit FileLock(const bf::path& path);
  /// Delete the lock file when the guard is destroyed.
  ~FileLock();
  /** Try to acquire the file lock.
   *
   *  \return Returns false if the file already exists.
   *          Returns true and creates the file otherwise.
   *
   *  \throws std::runtime_error if called another time after returning
   *  true or if the lockfile cannot be closed.
   */
  bool acquire();

 private:
  /// Path to the file lock.
  bf::path path_;
  /// Whether the lock has been acquired.
  bool acquired_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_FILELOCK_H_
