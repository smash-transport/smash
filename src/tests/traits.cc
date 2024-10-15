/*
 *
 *    Copyright (c) 2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/traits.h"

#include <array>
#include <map>
#include <optional>
#include <string>
#include <vector>

using namespace smash;

TEST(streamable) {
  std::vector<bool> v{
      is_writable_to_stream<std::stringstream, int>::value,
      is_writable_to_stream<std::stringstream, float>::value,
      is_writable_to_stream<std::stringstream, double>::value,
      is_writable_to_stream<std::stringstream, std::string>::value};
  for (const bool b : v)
    VERIFY(b);
}

TEST(not_streamable) {
  std::vector<bool> v{
      is_writable_to_stream<std::stringstream, std::vector<int>>::value,
      is_writable_to_stream<std::stringstream, std::array<double, 3>>::value,
      is_writable_to_stream<std::stringstream, std::optional<int>>::value,
      is_writable_to_stream<std::stringstream, std::map<int, double>>::value};
  for (const bool b : v)
    VERIFY(!b);
}
