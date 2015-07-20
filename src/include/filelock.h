/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include <boost/filesystem.hpp>

#include "forwarddeclarations.h"

namespace Smash {

/// Guard to create a file lock.
///
/// This can be used to make sure concurrent processes do not interfere when for
/// example writing files.
///
/// The lock file will be deleted when the FileLock is destroyed.
///
/// Internally this uses some POSIX system calls to atomically check for the
/// existence of a file and create it if it does not exists. This works even
/// over NFS (for at least NFSv3 and Linux 2.6).
class FileLock {
 public:
  /// Construct a file lock guard with a lock file at the given path.
  ///
  /// This will not create a file, use acquire() for that.
  FileLock(bf::path path);
  /// Delete the lock file when the guard is destroyed.
  ~FileLock();
  /// Try to acquire the file lock.
  ///
  /// Returns false if the file already exists.
  /// Returns true and creates the file otherwise.
  ///
  /// Will throw an std::runtime_error if called another time after returning
  /// true or if the lockfile cannot be closed.
  bool acquire();
  bf::path path_;

 private:
  bool acquired_;
};

}  // namespace Smash
