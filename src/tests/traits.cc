/*
 *
 *    Copyright (c) 2024-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/traits.h"

#include <array>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "smash/forwarddeclarations.h"

using namespace smash;

namespace std {

// Make set being streamable for this test file
template <typename T>
static std::ostream &operator<<(std::ostream &s, const std::set<T> &set) {
  s << '{';
  for (const auto &x : set) {
    s << x << ' ';
  }
  return s << '}';
}
// Make map being streamable for this test file but not pair to demonstrate
// that our trait implementation consider map writable to stream if its two
// types are, even if pair is not streamable.
template <typename K, typename V>
static std::ostream &operator<<(std::ostream &s, const std::map<K, V> &map) {
  s << '{';
  for (const auto &[key, value] : map) {
    s << '(' << key << ": " << value << ')';
  }
  return s << '}';
}

}  // namespace std

TEST(stl_container) {
  std::vector<bool> fulfilling{
      is_stl_container_v<std::vector<int>>, is_stl_container_v<std::set<float>>,
      is_stl_container_v<std::map<std::string, double>>,
      is_stl_container_v<std::map<ThermodynamicQuantity, double>>};
  for (const bool yes : fulfilling) {
    VERIFY(yes);
  }
  std::vector<bool> not_fulfilling{
      is_stl_container_v<int>, is_stl_container_v<ThermodynamicQuantity>,
      is_stl_container_v<std::array<int, 3>>,
      is_stl_container_v<std::unique_ptr<std::vector<int>>>};
  for (const bool no : not_fulfilling) {
    VERIFY(!no);
  }
}

TEST(tuple_like) {
  std::vector<bool> fulfilling{
      is_tuple_like_v<std::pair<int, int>>,
      is_tuple_like_v<std::array<int, 42>>,
      is_tuple_like_v<std::tuple<int, float, ThermodynamicQuantity>>};
  for (const bool yes : fulfilling) {
    VERIFY(yes);
  }
  std::vector<bool> not_fulfilling{
      is_tuple_like_v<int>, is_tuple_like_v<ThermodynamicQuantity>,
      is_tuple_like_v<std::vector<int>>,
      is_tuple_like_v<std::unique_ptr<std::pair<int, int>>>};
  for (const bool no : not_fulfilling) {
    VERIFY(!no);
  }
}

TEST(streamable) {
  std::vector<bool> fulfilling{
      is_streamable_v<std::stringstream, int>,
      is_streamable_v<std::stringstream, float>,
      is_streamable_v<std::stringstream, double>,
      is_streamable_v<std::stringstream, std::set<double>>,
      is_streamable_v<std::stringstream, std::string>,
      is_streamable_v<std::stringstream, std::map<int, double>>};
  for (const bool yes : fulfilling) {
    VERIFY(yes);
  }
  std::vector<bool> not_fulfilling{
      is_streamable_v<std::ostream, ThermodynamicQuantity>,
      is_streamable_v<std::ostream, std::vector<int>>,
      is_streamable_v<std::ostream, std::array<double, 3>>,
      is_streamable_v<std::ostream, std::optional<int>>};
  for (const bool no : not_fulfilling) {
    VERIFY(!no);
  }
}

TEST(writable_to_stream) {
  std::vector<bool> fulfilling{
      is_writable_to_stream_v<std::stringstream, int>,
      is_writable_to_stream_v<std::stringstream, std::string>,
      is_writable_to_stream_v<std::ostream, std::map<int, double>>,
      is_writable_to_stream_v<std::stringstream, std::set<double>>,
      is_writable_to_stream_v<std::stringstream, std::set<std::set<double>>>};
  for (const bool yes : fulfilling) {
    VERIFY(yes);
  }
  std::vector<bool> not_fulfilling{
      is_writable_to_stream_v<std::ostream, ThermodynamicQuantity>,
      is_writable_to_stream_v<std::ostream, std::vector<int>>,
      is_writable_to_stream_v<std::ostream, std::array<double, 3>>,
      is_writable_to_stream_v<std::ostream, std::optional<int>>,
      is_writable_to_stream_v<std::ostream, std::set<ThermodynamicQuantity>>};
  for (const bool no : not_fulfilling) {
    VERIFY(!no);
  }
}
