#include "include/smash/sha256.h"
#include <cstring>

// Macros

#define ror(value, bits) (((value) >> (bits)) | ((value) << (32 - (bits))))

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define STORE32H(x, y)                     \
  {                                        \
    (y)[0] = (uint8_t)(((x) >> 24) & 255); \
    (y)[1] = (uint8_t)(((x) >> 16) & 255); \
    (y)[2] = (uint8_t)(((x) >> 8) & 255);  \
    (y)[3] = (uint8_t)((x)&255);           \
  }

#define LOAD32H(x, y)                                                         \
  {                                                                           \
    x = ((uint32_t)((y)[0] & 255) << 24) | ((uint32_t)((y)[1] & 255) << 16) | \
        ((uint32_t)((y)[2] & 255) << 8) | ((uint32_t)((y)[3] & 255));         \
  }

#define STORE64H(x, y)                     \
  {                                        \
    (y)[0] = (uint8_t)(((x) >> 56) & 255); \
    (y)[1] = (uint8_t)(((x) >> 48) & 255); \
    (y)[2] = (uint8_t)(((x) >> 40) & 255); \
    (y)[3] = (uint8_t)(((x) >> 32) & 255); \
    (y)[4] = (uint8_t)(((x) >> 24) & 255); \
    (y)[5] = (uint8_t)(((x) >> 16) & 255); \
    (y)[6] = (uint8_t)(((x) >> 8) & 255);  \
    (y)[7] = (uint8_t)((x)&255);           \
  }

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

#define BLOCK_SIZE 64

// Internal functions

// Various logical functions
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

/**
 *  Compress 512-bits
 */
static void transform_function(Sha256Context* context, uint8_t const* buffer) {
  uint32_t S[8];
  uint32_t W[64];
  uint32_t t0;
  uint32_t t1;
  uint32_t t;
  int i;

  // Copy state into S
  for (i = 0; i < 8; i++) {
    S[i] = context->state[i];
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
    context->state[i] = context->state[i] + S[i];
  }
}

// Public functions

/// Initialises a SHA256 Context. Use this to initialise/reset a context.
void sha256_initialize(Sha256Context* context) {
  context->curlen = 0;
  context->length = 0;
  context->state[0] = 0x6A09E667UL;
  context->state[1] = 0xBB67AE85UL;
  context->state[2] = 0x3C6EF372UL;
  context->state[3] = 0xA54FF53AUL;
  context->state[4] = 0x510E527FUL;
  context->state[5] = 0x9B05688CUL;
  context->state[6] = 0x1F83D9ABUL;
  context->state[7] = 0x5BE0CD19UL;
}

/**
 * Adds data to the SHA256 context. This will process the data and update the
 * internal state of the context. Keep on calling this function until all the
 * data has been added. Then call sha256_finalize to calculate the hash.
 */
void sha256_update(Sha256Context* context, uint8_t const* buffer,
                  uint32_t buffer_size) {
  uint32_t n;

  if (context->curlen > sizeof(context->buf)) {
    return;
  }

  while (buffer_size > 0) {
    if (context->curlen == 0 && buffer_size >= BLOCK_SIZE) {
      transform_function(context, (uint8_t*)buffer);
      context->length += BLOCK_SIZE * 8;
      buffer = (uint8_t*)buffer + BLOCK_SIZE;
      buffer_size -= BLOCK_SIZE;
    } else {
      n = MIN(buffer_size, (BLOCK_SIZE - context->curlen));
      memcpy(context->buf + context->curlen, buffer, (size_t)n);
      context->curlen += n;
      buffer = (uint8_t*)buffer + n;
      buffer_size -= n;
      if (context->curlen == BLOCK_SIZE) {
        transform_function(context, context->buf);
        context->length += 8 * BLOCK_SIZE;
        context->curlen = 0;
      }
    }
  }
}

/**
 *  Performs the final calculation of the hash and returns the digest (32 byte
 * buffer containing 256bit hash). After calling this, sha256_initialized must be
 * used to reuse the context.
 */
void sha256_finalize(Sha256Context* context, SHA256_HASH* digest) {
  int i;

  if (context->curlen >= sizeof(context->buf)) {
    return;
  }

  // Increase the length of the message
  context->length += context->curlen * 8;

  // Append the '1' bit
  context->buf[context->curlen++] = (uint8_t)0x80;

  // if the length is currently above 56 bytes we append zeros
  // then compress.  Then we can fall back to padding zeros and length
  // encoding like normal.
  if (context->curlen > 56) {
    while (context->curlen < 64) {
      context->buf[context->curlen++] = (uint8_t)0;
    }
    transform_function(context, context->buf);
    context->curlen = 0;
  }

  // Pad up to 56 bytes of zeroes
  while (context->curlen < 56) {
    context->buf[context->curlen++] = (uint8_t)0;
  }

  // Store length
  STORE64H(context->length, context->buf + 56);
  transform_function(context, context->buf);

  // Copy output
  for (i = 0; i < 8; i++) {
    STORE32H(context->state[i], digest->data() + (4 * i));
  }
}

/**
 *  Combines sha256_initialize, sha256_update, and sha256_finalize into one
 * function. Calculates the SHA256 hash of the buffer.
 */
void sha256_calculate(uint8_t const* buffer, uint32_t buffer_size,
                      SHA256_HASH* digest) {
  Sha256Context context;

  sha256_initialize(&context);
  sha256_update(&context, buffer, buffer_size);
  sha256_finalize(&context, digest);
}
