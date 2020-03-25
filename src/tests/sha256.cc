/*
 *
 *    Copyright (c) 2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <vir/test.h>  // This include has to be first

#include <string>

#include "../include/smash/sha256.h"

using namespace smash;

struct TestVector {
  std::vector<uint8_t> input;
  std::string output;
};

static const std::vector<TestVector> test_vectors = {
    {{0xbd},
     "68325720aabd7c82f30f554b313d0570c95accbb7dc4b5aae11204c08ffe732b"},
    {{0xc9, 0x8c, 0x8e, 0x55},
     "7abc22c0ae5af26ce93dbb94433a0e0b2e119d014f8e7f65bd56c61ccccd9504"}};

TEST(basic) {
  for (const auto& t : test_vectors) {
    const auto digest = sha256::calculate(t.input.data(), t.input.size());
    COMPARE(t.output, sha256::hash_to_string(digest));
  }
}
