#pragma once

#include <cstdint>
#include <cstdio>
#include <array>

struct Sha256Context {
  uint64_t length;
  uint32_t state[8];
  uint32_t curlen;
  uint8_t buf[64];
};

#define SHA256_HASH_SIZE (256 / 8)

typedef std::array<uint8_t, SHA256_HASH_SIZE> SHA256_HASH;

// Public functions

/// Initialises a SHA256 Context. Use this to initialise/reset a context.
void sha256_initialise(Sha256Context* context);

/**
 * Adds data to the SHA256 context. This will process the data and update the
 * internal state of the context. Keep on calling this function until all the
 * data has been added. Then call sha256_finalize to calculate the hash.
 */
void sha256_update(Sha256Context* context, void const* buffer,
                   uint32_t buffer_size);

/**
 * Performs the final calculation of the hash and returns the digest (32 byte
 * buffer containing 256bit hash). After calling this, sha256_initialized must
 * be used to reuse the context.
 */
void sha256_finalize(Sha256Context* context, SHA256_HASH* digest);

/**
 *  Combines sha256_initialize, sha256_update, and sha256_finalize into one
 * function. Calculates the SHA256 hash of the buffer.
 */
void sha256_calculate(uint8_t const* buffer, uint32_t buffer_size,
                      SHA256_HASH* digest);
