/*
 *    Copyright (c) 2019-2020,2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

// This is based on the public-domain implementation from
// https://github.com/WaterJuice/WjCryptLib.

#include "smash/sha256.h"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>

// Macros

#define ror(value, bits) (((value) >> (bits)) | ((value) << (32 - (bits))))

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define STORE32H(x, y)                                \
  {                                                   \
    (y)[0] = static_cast<uint8_t>(((x) >> 24) & 255); \
    (y)[1] = static_cast<uint8_t>(((x) >> 16) & 255); \
    (y)[2] = static_cast<uint8_t>(((x) >> 8) & 255);  \
    (y)[3] = static_cast<uint8_t>((x)&255);           \
  }

#define LOAD32H(x, y)                                 \
  {                                                   \
    x = (static_cast<uint32_t>((y)[0] & 255) << 24) | \
        (static_cast<uint32_t>((y)[1] & 255) << 16) | \
        (static_cast<uint32_t>((y)[2] & 255) << 8) |  \
        (static_cast<uint32_t>((y)[3] & 255));        \
  }

#define STORE64H(x, y)                                \
  {                                                   \
    (y)[0] = static_cast<uint8_t>(((x) >> 56) & 255); \
    (y)[1] = static_cast<uint8_t>(((x) >> 48) & 255); \
    (y)[2] = static_cast<uint8_t>(((x) >> 40) & 255); \
    (y)[3] = static_cast<uint8_t>(((x) >> 32) & 255); \
    (y)[4] = static_cast<uint8_t>(((x) >> 24) & 255); \
    (y)[5] = static_cast<uint8_t>(((x) >> 16) & 255); \
    (y)[6] = static_cast<uint8_t>(((x) >> 8) & 255);  \
    (y)[7] = static_cast<uint8_t>((x)&255);           \
  }

#define Ch(x, y, z) (z ^ (x & (y ^ z)))
#define Maj(x, y, z) (((x | y) & z) | (x & y))
#define S(x, n) ror((x), (n))
#define R(x, n) (((x)&0xFFFFFFFFUL) >> (n))
#define Sigma0(x) (S(x, 2) ^ S(x, 13) ^ S(x, 22))
#define Sigma1(x) (S(x, 6) ^ S(x, 11) ^ S(x, 25))
#define Gamma0(x) (S(x, 7) ^ S(x, 18) ^ R(x, 3))
#define Gamma1(x) (S(x, 17) ^ S(x, 19) ^ R(x, 10))

#define Sha256Round(a, b, c, d, e, f, g, h, i)    \
  t0 = h + Sigma1(e) + Ch(e, f, g) + K[i] + W[i]; \
  t1 = Sigma0(a) + Maj(a, b, c);                  \
  d += t0;                                        \
  h = t0 + t1;

// Constants

/// The K array.
static const uint32_t K[64] = {
    0x428a2f98UL, 0x71374491UL, 0xb5c0fbcfUL, 0xe9b5dba5UL, 0x3956c25bUL,
    0x59f111f1UL, 0x923f82a4UL, 0xab1c5ed5UL, 0xd807aa98UL, 0x12835b01UL,
    0x243185beUL, 0x550c7dc3UL, 0x72be5d74UL, 0x80deb1feUL, 0x9bdc06a7UL,
    0xc19bf174UL, 0xe49b69c1UL, 0xefbe4786UL, 0x0fc19dc6UL, 0x240ca1ccUL,
    0x2de92c6fUL, 0x4a7484aaUL, 0x5cb0a9dcUL, 0x76f988daUL, 0x983e5152UL,
    0xa831c66dUL, 0xb00327c8UL, 0xbf597fc7UL, 0xc6e00bf3UL, 0xd5a79147UL,
    0x06ca6351UL, 0x14292967UL, 0x27b70a85UL, 0x2e1b2138UL, 0x4d2c6dfcUL,
    0x53380d13UL, 0x650a7354UL, 0x766a0abbUL, 0x81c2c92eUL, 0x92722c85UL,
    0xa2bfe8a1UL, 0xa81a664bUL, 0xc24b8b70UL, 0xc76c51a3UL, 0xd192e819UL,
    0xd6990624UL, 0xf40e3585UL, 0x106aa070UL, 0x19a4c116UL, 0x1e376c08UL,
    0x2748774cUL, 0x34b0bcb5UL, 0x391c0cb3UL, 0x4ed8aa4aUL, 0x5b9cca4fUL,
    0x682e6ff3UL, 0x748f82eeUL, 0x78a5636fUL, 0x84c87814UL, 0x8cc70208UL,
    0x90befffaUL, 0xa4506cebUL, 0xbef9a3f7UL, 0xc67178f2UL};

namespace smash {
namespace sha256 {

constexpr size_t BLOCK_SIZE = 64;

// Internal functions

/**
 *  Compress 512-bits
 */
void Context::transform_function(uint8_t const* buffer) {
  uint32_t S[8];
  uint32_t W[64];
  uint32_t t0;
  uint32_t t1;
  uint32_t t;
  int i;

  // Copy state into S
  for (i = 0; i < 8; i++) {
    S[i] = state_[i];
  }

  // Copy the state into 512-bits into W[0..15]
  for (i = 0; i < 16; i++) {
    LOAD32H(W[i], buffer + (4 * i));
  }

  // Fill W[16..63]
  for (i = 16; i < 64; i++) {
    W[i] = Gamma1(W[i - 2]) + W[i - 7] + Gamma0(W[i - 15]) + W[i - 16];
  }

  // Compress
  for (i = 0; i < 64; i++) {
    Sha256Round(S[0], S[1], S[2], S[3], S[4], S[5], S[6], S[7], i);
    t = S[7];
    S[7] = S[6];
    S[6] = S[5];
    S[5] = S[4];
    S[4] = S[3];
    S[3] = S[2];
    S[2] = S[1];
    S[1] = S[0];
    S[0] = t;
  }

  // Feedback
  for (i = 0; i < 8; i++) {
    state_[i] += S[i];
  }
}

// Public functions

void Context::reset() {
  curlen_ = 0;
  length_ = 0;
  state_[0] = 0x6A09E667UL;
  state_[1] = 0xBB67AE85UL;
  state_[2] = 0x3C6EF372UL;
  state_[3] = 0xA54FF53AUL;
  state_[4] = 0x510E527FUL;
  state_[5] = 0x9B05688CUL;
  state_[6] = 0x1F83D9ABUL;
  state_[7] = 0x5BE0CD19UL;
}

void Context::update(uint8_t const* buffer, size_t buffer_size) {
  size_t n;

  if (curlen_ > sizeof(buf_)) {
    return;
  }

  while (buffer_size > 0) {
    if (curlen_ == 0 && buffer_size >= BLOCK_SIZE) {
      transform_function(buffer);
      length_ += BLOCK_SIZE * 8;
      buffer += BLOCK_SIZE;
      buffer_size -= BLOCK_SIZE;
    } else {
      n = MIN(buffer_size, (BLOCK_SIZE - curlen_));
      std::memcpy(buf_ + curlen_, buffer, n);
      curlen_ += n;
      buffer += n;
      buffer_size -= n;
      if (curlen_ == BLOCK_SIZE) {
        transform_function(buf_);
        length_ += 8 * BLOCK_SIZE;
        curlen_ = 0;
      }
    }
  }
}

void Context::update(const std::string& buffer) {
  update(reinterpret_cast<const uint8_t*>(buffer.c_str()), buffer.size());
}

Hash Context::finalize() {
  Hash digest{{}};
  if (curlen_ >= sizeof(buf_)) {
    return digest;
  }

  // Increase the length of the message
  length_ += curlen_ * 8;

  // Append the '1' bit
  buf_[curlen_++] = 0x80;

  // if the length is currently above 56 bytes we append zeros
  // then compress.  Then we can fall back to padding zeros and length
  // encoding like normal.
  if (curlen_ > 56) {
    while (curlen_ < 64) {
      buf_[curlen_++] = 0;
    }
    transform_function(buf_);
    curlen_ = 0;
  }

  // Pad up to 56 bytes of zeroes
  while (curlen_ < 56) {
    buf_[curlen_++] = 0;
  }

  // Store length
  STORE64H(length_, buf_ + 56);
  transform_function(buf_);

  // Copy output
  for (int i = 0; i < 8; i++) {
    STORE32H(state_[i], digest.data() + (4 * i));
  }
  return digest;
}

Hash calculate(uint8_t const* buffer, size_t buffer_size) {
  Context context;
  context.update(buffer, buffer_size);
  return context.finalize();
}

std::string hash_to_string(Hash hash) {
  std::stringstream ss;
  ss << std::hex;
  for (uint16_t i : hash) {
    ss << std::setw(2) << std::setfill('0') << i;
  }
  return ss.str();
}

}  // namespace sha256
}  // namespace smash
