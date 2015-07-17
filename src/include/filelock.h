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

class FileLock {
 public:
  FileLock(bf::path path);
  ~FileLock();
  bool acquire();
  bf::path path_;

 private:
  bool acquired_;
};

}  // namespace Smash
