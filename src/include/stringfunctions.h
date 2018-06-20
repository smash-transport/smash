/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_STRINGFUNCTIONS_H_
#define SRC_INCLUDE_STRINGFUNCTIONS_H_

#include <string>
#include <vector>

namespace smash {

/**
 * Strip leading and trailing whitespaces.
 *
 * \param s String to be trimmed.
 * \return Trimmed string.
 */
std::string trim(const std::string &s);

/**
 * Remove all instances of a substring p in a string s.
 *
 * \param[inout] s String to be searched and modified.
 * \param[in] p Substring to be removed.
 */
void remove_substr(std::string &s, const std::string &p);

/**
 * Remove ⁺, ⁻, ⁰ from string.
 * 
 * \param[inout] s String to be cleaned.
 */
void isoclean(std::string &s);

/**
 * Split string by delimiter.
 *
 * \param[in] s String to be split.
 * \param[in] delim Splitting delimiter.
 * \return Split string.
 */
std::vector<std::string> split(const std::string &s, char delim);

namespace utf8 {
/**
 * Fill string with characters to the left until the given width is reached.
 *
 * \param[in] s Input string.
 * \param[in] width Total width of output string.
 * \param[in] fill Filling character.
 * \return Padded string.
 */
std::string fill_left(const std::string &s, size_t width, char fill = ' ');

/**
 * Fill string with characters to the right until the given width is reached.
 *
 * \param[in] s Input string.
 * \param[in] width Total width of output string.
 * \param[in] fill Filling character.
 * \return Padded string.
 */
std::string fill_right(const std::string &s, size_t width, char fill = ' ');

/**
 * Fill string with characters at both sides until the given width is reached.
 *
 * \param[in] s Input string.
 * \param[in] width Total width of output string.
 * \param[in] fill Filling character.
 * \return Padded string.
 */
std::string fill_both(const std::string &s, size_t width, char fill = ' ');

/**
 * Extract the first byte from a given value.
 *
 * This function was taken from the Boost-licensed library UTF8-CPP.
 * See http://utfcpp.sourceforge.net/.
 *
 * \tparam octet_type Type for one byte
 */
template <typename octet_type>
inline uint8_t mask8(octet_type oc) {
  return static_cast<uint8_t>(0xff & oc);
}

/**
 * Given an iterator to the beginning of a UTF-8 sequence, return the
 * length of the next UTF-8 code point.
 *
 * This function was taken from the Boost-licensed library UTF8-CPP.
 * See http://utfcpp.sourceforge.net/.
 */
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

}  // namespace smash

#endif  // SRC_INCLUDE_STRINGFUNCTIONS_H_
