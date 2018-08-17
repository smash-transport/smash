#include "unittest.h"  // This include has to be first

#include "../include/smash/stringfunctions.h"

using namespace smash::utf8;

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
