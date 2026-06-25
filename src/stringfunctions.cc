/*
 *
 *    Copyright (c) 2014-2018,2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/stringfunctions.h"

#include <numeric>
#include <sstream>

namespace smash {

namespace utf8 {

/**
 * Adjust filling width by taking the size of unicode characters into account.
 * This is necessary, because UTF-8 characters can be represented by more than
 * byte.
 *
 * \param[in] s String to be filled.
 * \param[in] width Width (in bytes) to be adjusted.
 * \return Adjusted width.
 */
inline static size_t adjust(const std::string &s, size_t width) {
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

std::string fill_left(const std::string &s, size_t width, char fill) {
  width = adjust(s, width - s.size());
  if (width > 0) {
    return std::string(width, fill) + s;
  }
  return s;
}

std::string fill_right(const std::string &s, size_t width, char fill) {
  width = adjust(s, width - s.size());
  if (width > 0) {
    return s + std::string(width, fill);
  }
  return s;
}

std::string fill_both(const std::string &s, size_t width, char fill) {
  width = adjust(s, width - s.size());
  if (width > 0) {
    const int l = width / 2;
    const int r = width - l;
    return std::string(l, fill) + s + std::string(r, fill);
  }
  return s;
}

}  // namespace utf8

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

/**
 * Split string by delimiter.
 *
 * \param[in] s String to be split.
 * \param[in] delim Splitting delimiter.
 * \param[out] result Split string as iterator.
 *
 * Necessary for the next function
 */
template <typename Out>
void split(const std::string &s, char delim, Out result);

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

template <typename Container>
static std::string join_impl(const Container &container,
                             std::string_view delim) {
  // Enable ADL (Argument-Dependent Lookup); not necessary now but harmless
  using std::begin;
  using std::end;
  using value_type = std::decay_t<decltype(*begin(container))>;
  /* Here this is not really needed as the developer calls this function with
   * the appropriate types, but it is ready to be moved to the header one day */
  static_assert(
      std::is_convertible_v<value_type, std::string_view>,
      "join() requires string-like objects convertible to std::string_view!");
  auto it = begin(container);
  const auto last = end(container);
  if (it == last) {
    return {};
  }
  std::size_t total_size = 0;
  std::size_t count = 0;
  for (const auto &s : container) {
    total_size += std::string_view{s}.size();
    ++count;
  }
  total_size += delim.size() * (count - 1);
  std::string result{};
  result.reserve(total_size);
  result.append(std::string_view{*it});
  ++it;
  for (; it != last; ++it) {
    result.append(delim);
    result.append(std::string_view{*it});
  }
  return result;
}

std::string join(const std::vector<std::string> &v, std::string_view delim) {
  return join_impl(v, delim);
}

std::string join(const std::vector<std::string_view> &v,
                 std::string_view delim) {
  return join_impl(v, delim);
}

std::string join(const std::set<std::string> &s, std::string_view delim) {
  return join_impl(s, delim);
}

std::string quote(const std::string &s) { return "\"" + s + "\""; }

}  // namespace smash
