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

// std::string fill_left(const std::string &s, int width, char fill = ' ');
std::string fill_right(const std::string &s, int width, char fill = ' ');
std::string fill_both(const std::string &s, int width, char fill = ' ');

/// takes a string and strips leading and trailing whitespaces.
std::string trim(const std::string &s);

namespace utf8 {
    // The functions here were taken from the Boost-licensed library UTF8-CPP.
    // See http://utfcpp.sourceforge.net/.

    /// Extract the first byte from a given value.
    template <typename octet_type>
    inline uint8_t mask8(octet_type oc) {
      return static_cast<uint8_t>(0xff & oc);
    }

    /// Given an iterator to the beginning of a UTF-8 sequence, return the
    /// length of the next UTF-8 code point.
    template <typename octet_iterator>
    inline typename std::iterator_traits<octet_iterator>::difference_type
    sequence_length(octet_iterator lead_it) {
      uint8_t lead = mask8(*lead_it);
      if (lead < 0x80)
        return 1;
      else if ((lead >> 5) == 0x6)
        return 2;
      else if ((lead >> 4) == 0xe)
        return 3;
      else if ((lead >> 3) == 0x1e)
        return 4;
      else
        return 0;
    }

}  // namespace utf8

}  // namespace Smash

#endif  // SRC_INCLUDE_STRINGFUNCTIONS_H_
