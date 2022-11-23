/*
 *
 *    Copyright (c) 2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/file.h"

namespace smash {

FilePtr fopen(const std::filesystem::path& filename, const std::string& mode) {
  FilePtr f{std::fopen(filename.c_str(), mode.c_str())};
  return f;
}

RenamingFilePtr::RenamingFilePtr(const std::filesystem::path& filename,
                                 const std::string& mode) {
  filename_ = filename;
  filename_unfinished_ = filename;
  filename_unfinished_ += ".unfinished";
  file_ = std::fopen(filename_unfinished_.c_str(), mode.c_str());
}

FILE* RenamingFilePtr::get() { return file_; }

RenamingFilePtr::~RenamingFilePtr() {
  std::fclose(file_);
  // we rename the output file only if we are not unwinding the stack
  // because of an exception
  if (std::uncaught_exceptions() == 0) {
    std::filesystem::rename(filename_unfinished_, filename_);
  }
}

}  // namespace smash
