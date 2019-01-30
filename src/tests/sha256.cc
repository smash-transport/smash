/*
 *
 *    Copyright (c) 2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include <string>

#include "../include/smash/sha256.h"

using namespace smash;

static std::string hash_to_string(sha256::Hash hash) {
  std::stringstream ss;
  ss << std::hex;
  for (uint16_t i : hash) {
    ss << std::setw(2) << std::setfill('0') << i;
  }
  return ss.str();
}

struct TestVector {
  std::vector<uint8_t> input;
  std::string output;
};

static const std::vector<TestVector> test_vectors = {
  { { 0xbd }, "68325720aabd7c82f30f554b313d0570c95accbb7dc4b5aae11204c08ffe732b" },
  { { 0xc9, 0x8c, 0x8e, 0x55 }, "7abc22c0ae5af26ce93dbb94433a0e0b2e119d014f8e7f65bd56c61ccccd9504" }
};

TEST(basic) {
  for (const auto& t : test_vectors) {
    sha256::Hash digest;
    sha256::calculate(t.input.data(), t.input.size(), &digest);
    COMPARE(t.output, hash_to_string(digest));
  }
}
