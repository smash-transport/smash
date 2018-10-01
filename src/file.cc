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

FilePtr fopen(const bf::path& filename, const std::string& mode) {
  FilePtr f{std::fopen(filename.c_str(), mode.c_str())};
  return f;
}

RenamingFilePtr::RenamingFilePtr(const bf::path& filename,
                                 const std::string& mode) {
  filename_ = filename;
  filename_unfinished_ = filename;
  filename_unfinished_ += ".unfinished";
  file_ = std::fopen(filename_unfinished_.c_str(), mode.c_str());
}

FILE* RenamingFilePtr::get() { return file_; }

RenamingFilePtr::~RenamingFilePtr() {
  std::fclose(file_);
  bf::rename(filename_unfinished_, filename_);
}

}  // namespace smash
