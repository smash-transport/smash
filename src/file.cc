/*
 *
 *    Copyright (c) 2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/file.h"

namespace smash {

FilePtr fopen(const std::string& filename, const std::string& mode) {
    FilePtr f{ std::fopen(filename.c_str(), mode.c_str()) };
    return f;
}

}  // namespace smash
