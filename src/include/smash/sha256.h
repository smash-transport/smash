#ifndef SRC_INCLUDE_SHA256_H_
#define SRC_INCLUDE_SHA256_H_

#include <cstdint>
#include <cstdio>
#include <array>

namespace smash {
namespace sha256 {

constexpr size_t HASH_SIZE = 256 / 8;

struct Context {
  uint64_t length;
  uint32_t state[8];
  uint32_t curlen;
  uint8_t buf[64];
};


typedef std::array<uint8_t, HASH_SIZE> Hash;

// Public functions

/// Initializes a SHA256 Context. Use this to initialize/reset a context.
void initialize(Context* context);

/**
 * Adds data to the SHA256 context. This will process the data and update the
 * internal state of the context. Keep on calling this function until all the
 * data has been added. Then call finalize to calculate the hash.
 */
void update(Context* context, uint8_t const* buffer,
            uint32_t buffer_size);

/**
 * Performs the final calculation of the hash and returns the digest (32 byte
 * buffer containing 256bit hash). After calling this, initialize must
 * be used to reuse the context.
 */
void finalize(Context* context, Hash* digest);

/**
 * Calculate the SHA256 hash of the buffer.
 */
void calculate(uint8_t const* buffer, uint32_t buffer_size,
               Hash* digest);

}  // namespace sha256
}  // namespace smash

#endif  // SRC_INCLUDE_SHA256_H_
