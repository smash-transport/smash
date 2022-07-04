/*
 *
 *    Copyright (c) 2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/numeric_cast.h"

TEST_CATCH(expect_overflow, std::overflow_error) {
  const size_t input = -1;
  smash::numeric_cast<uint32_t>(input);
}

TEST(successful_conversion) {
  const size_t input = 12345;
  const auto output = smash::numeric_cast<uint32_t>(input);
  COMPARE(input, output);
}