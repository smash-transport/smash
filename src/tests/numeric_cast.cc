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

#include <limits>

TEST_CATCH(positive_overflow_signed_to_signed, std::domain_error) {
  const int input = std::numeric_limits<int16_t>::max() + 1;
  smash::numeric_cast<int16_t>(input);
}
TEST_CATCH(positive_overflow_signed_to_unsigned, std::domain_error) {
  const int input = std::numeric_limits<uint16_t>::max() + 1;
  smash::numeric_cast<uint16_t>(input);
}
TEST_CATCH(positive_overflow_unsigned_to_signed, std::domain_error) {
  const unsigned int input = std::numeric_limits<uint16_t>::max() + 1;
  smash::numeric_cast<int16_t>(input);
}
TEST_CATCH(positive_overflow_unsigned_to_unsigned, std::domain_error) {
  const unsigned int input = std::numeric_limits<uint16_t>::max() + 1;
  smash::numeric_cast<uint16_t>(input);
}

TEST_CATCH(negative_overflow_signed_to_signed, std::domain_error) {
  const int input = std::numeric_limits<int16_t>::lowest() - 1;
  smash::numeric_cast<int16_t>(input);
}
TEST_CATCH(negative_overflow_signed_to_unsigned, std::domain_error) {
  const size_t input = -1;
  smash::numeric_cast<uint32_t>(input);
}

TEST_CATCH(non_representable_double_to_double, std::domain_error) {
  const uint64_t input = std::numeric_limits<int64_t>::max();
  smash::numeric_cast<double>(input);
}

TEST_CATCH(double_to_integer_loosing_conversion, std::domain_error) {
  const double input = -3.14;
  smash::numeric_cast<int>(input);
}

TEST(successful_conversion) {
  const size_t input = 12345;
  const auto output = smash::numeric_cast<uint32_t>(input);
  COMPARE(input, output);
}

TEST(successful_dangerous_conversion) {
  const double input = 42.0;
  const auto output = smash::numeric_cast<int>(input);
  COMPARE(input, output);
}
