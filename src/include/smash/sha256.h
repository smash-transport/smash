#ifndef SRC_INCLUDE_SHA256_H_
#define SRC_INCLUDE_SHA256_H_

#include <array>
#include <cstdint>
#include <cstdio>
#include <string>

namespace smash {
namespace sha256 {

/// Size of a SHA256 hash.
constexpr size_t HASH_SIZE = 256 / 8;

/// A SHA256 hash.
typedef std::array<uint8_t, HASH_SIZE> Hash;

/// A SHA256 context.
class Context {
 private:
  /// Length of the SHA256 hash.
  uint64_t length_;
  /// State of the SHA256 hash.
  uint32_t state_[8];
  /// Current length of the SHA256 hash.
  size_t curlen_;
  /// Buffer of the SHA256 hash.
  uint8_t buf_[64];

  /**
   *  Compress 512-bits
   */
  void transform_function(uint8_t const* buffer);

 public:
  /// Construct a SHA256 context.
  Context() { reset(); }

  /// Reset the SHA256 context.
  void reset();

  /**
   * Add data to the SHA256 context. This will process the data and update the
   * internal state of the context. Keep on calling this function until all the
   * data has been added. Then call `finalize` to calculate the hash.
   */
  void update(uint8_t const* buffer, size_t buffer_size);

  /**
   * Add data to the SHA256 context. This will process the data and update the
   * internal state of the context. Keep on calling this function until all the
   * data has been added. Then call `finalize` to calculate the hash.
   */
  void update(const std::string& buffer);

  /**
   * Performs the final calculation of the hash and returns the digest (32 byte
   * buffer containing 256bit hash). After calling this, `reset` must
   * be used to reuse the context.
   */
  Hash finalize();
};

/**
 * Calculate the SHA256 hash of the buffer.
 */
Hash calculate(uint8_t const* buffer, size_t buffer_size);

/**
 * Convert a SHA256 hash to a hexadecimal string.
 */
std::string hash_to_string(Hash hash);

}  // namespace sha256
}  // namespace smash

#endif  // SRC_INCLUDE_SHA256_H_
