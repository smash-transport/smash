/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_STRINGFUNCTIONS_H_
#define SRC_INCLUDE_STRINGFUNCTIONS_H_

#include <string>

namespace Smash {

std::string fill_left(const std::string &s, int width, char fill = ' ');
std::string fill_right(const std::string &s, int width, char fill = ' ');
std::string fill_both(const std::string &s, int width, char fill = ' ');

/// takes a string and strips leading and trailing whitespaces.
std::string trim(const std::string &s);

}  // namespace Smash

#endif  // SRC_INCLUDE_STRINGFUNCTIONS_H_
