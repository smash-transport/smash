/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/stringfunctions.h"

#include <sstream>

namespace smash {

inline static int utf8_adjust(const std::string &s, int width) {
  for (unsigned char c : s) {
    if (c >= 0xFC) {
      width += 5;
    } else if (c >= 0xF8) {
      width += 4;
    } else if (c >= 0xF0) {
      width += 3;
    } else if (c >= 0xE0) {
      width += 2;
    } else if (c == 0xCC || c == 0xCD) {
      // combining character (2 Bytes) - doesn't appear at all
      width += 2;
    } else if (c >= 0xC0) {
      width += 1;
    }
  }
  return width;
}

/*std::string fill_left(const std::string &s, int width, char fill) {
  width = utf8_adjust(s, width - s.size());
  if (width > 0) {
    return std::string(width, fill) + s;
  }
  return s;
}*/

std::string fill_right(const std::string &s, int width, char fill) {
  width = utf8_adjust(s, width - s.size());
  if (width > 0) {
    return s + std::string(width, fill);
  }
  return s;
}

std::string fill_both(const std::string &s, int width, char fill) {
  width = utf8_adjust(s, width - s.size());
  if (width > 0) {
    const int l = width / 2;
    const int r = width - l;
    return std::string(l, fill) + s + std::string(r, fill);
  }
  return s;
}

std::string trim(const std::string &s) {
  const auto begin = s.find_first_not_of(" \t\n\r");
  if (begin == std::string::npos) {
    return {};
  }
  const auto end = s.find_last_not_of(" \t\n\r");
  return s.substr(begin, end - begin + 1);
}

void remove_substr(std::string &s, const std::string &p) {
  using str = std::string;
  str::size_type n = p.length();
  for (str::size_type i = s.find(p); i != str::npos; i = s.find(p)) {
    s.erase(i, n);
  }
}

void isoclean(std::string &s) {
  remove_substr(s, "⁺");
  remove_substr(s, "⁻");
  remove_substr(s, "⁰");
}

template <typename Out>
void split(const std::string &s, char delim, Out result) {
  std::stringstream ss;
  ss.str(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    *(result++) = item;
  }
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, std::back_inserter(elems));
  return elems;
}

}  // namespace smash
