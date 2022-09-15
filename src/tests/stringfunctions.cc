/*
 *
 *    Copyright (c) 2015,2017-2020,2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/stringfunctions.h"

#include <sstream>
#include <vector>

namespace smash {
namespace utf8 {
TEST(sequence_length) {
  COMPARE(sequence_length("xπ"), 1);
  COMPARE(sequence_length("πx"), 2);
  COMPARE(sequence_length("ᛒ x"), 3);
  COMPARE(sequence_length("🅑 x"), 4);
}

TEST(fill) {
  COMPARE(fill_left("xπ", 5), "   xπ");
  COMPARE(fill_right("xπ", 5), "xπ   ");
  COMPARE(fill_both("xπ", 5), " xπ  ");

  COMPARE(fill_left("xπ", 5, '#'), "###xπ");
  COMPARE(fill_right("xπ", 5, '#'), "xπ###");
  COMPARE(fill_both("xπ", 5, '#'), "#xπ##");
}
}  // namespace utf8

TEST(trim_and_remove) {
  COMPARE(trim("  xπ# \n \r \t "), "xπ#");
  COMPARE(trim("  xπ# \n  xπ "), "xπ# \n  xπ");

  std::string s = "xπ#\n xπ90009_pion_\n";
  remove_substr(s, "xπ");
  COMPARE(s, "#\n 90009_pion_\n");
  remove_substr(s, "\n");
  COMPARE(s, "# 90009_pion_");
  remove_substr(s, " 90009_pion_");
  COMPARE(s, "#");
}

TEST(split) {
  std::string s = "xπ/.90928/_is_pion_/true";
  std::vector<std::string> storage;
  storage = split(s, '/');
  COMPARE(storage[0], "xπ");
  COMPARE(storage[1], ".90928");
  COMPARE(storage[2], "_is_pion_");
}

TEST(join) {
  const std::vector<std::string> v = {"Hello", "my", "SMASHies!"};
  const std::string j1 = "Hello my SMASHies!";
  const std::string j2 = "Hello__my__SMASHies!";
  COMPARE(join(v, " "), j1);
  COMPARE(join(v, "__"), j2);
}

TEST(quote) {
  const std::string input = "To be quoted";
  const std::string result = "\"To be quoted\"";
  COMPARE(quote(input), result);
}


}  // namespace smash
