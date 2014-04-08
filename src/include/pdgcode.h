/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_PDGCODE_H_
#define SRC_INCLUDE_PDGCODE_H_

#include <cstdio>
// #include <istringstream>
#include <stdexcept>
#include <string>

namespace Smash {

/** PDGCode stores a Particle Data Group Particle Numbering Scheme
 * particle type number.
 *
 * Usage:
 * ------
 * \code
 * #include "include/pdgcode.h"
 *
 * // needs to be filled.
 * \endcode
 *
 **/

class PDGCode {
 public:
  /// thrown for invalid values for theta
  struct InvalidPDGCode : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };

  /// Standard initializer
  //PDGCode() {}
  /** Initialize using a string
   *
   * The string is interpreted as a hexadecimal number, i.e., @211@ is
   * interpreted as @0x211 = 529_{10}@.
   */
  // PDGCode(std::stringstream codestring) {
  //   codestring >> std::hex >> pdgcode_;
  // }
  PDGCode(const int codenumber) {
    unsigned int bitmask = 0xf;
    int abscode = std::abs(codenumber);
    printf("bitmask: 0x%08x\nabscode: 0x%08x\nVALUE:   0x%08x\n", bitmask, abscode, abscode & bitmask);
    n_J_ = abscode & bitmask;
    abscode >>= 4;
    printf("bitmask: 0x%08x\nabscode: 0x%08x\nVALUE:   0x%08x\n", bitmask, abscode, abscode & bitmask);
    n_q3_ = abscode & bitmask;
    abscode >>= 4;
    printf("bitmask: 0x%08x\nabscode: 0x%08x\nVALUE:   0x%08x\n", bitmask, abscode, abscode & bitmask);
    n_q2_ = abscode & bitmask;
    abscode >>= 4;
    printf("bitmask: 0x%08x\nabscode: 0x%08x\nVALUE:   0x%08x\n", bitmask, abscode, abscode & bitmask);
    n_q1_ = abscode & bitmask;
    abscode >>= 4;
    printf("bitmask: 0x%08x\nabscode: 0x%08x\nVALUE:   0x%08x\n", bitmask, abscode, abscode & bitmask);
    n_L_  = abscode & bitmask;
    abscode >>= 4;
    printf("bitmask: 0x%08x\nabscode: 0x%08x\nVALUE:   0x%08x\n", bitmask, abscode, abscode & bitmask);
    n_R_  = abscode & bitmask;
    abscode >>= 4;
    printf("bitmask: 0x%08x\nabscode: 0x%08x\nVALUE:   0x%08x\n", bitmask, abscode, abscode & bitmask);
    n_    = abscode & bitmask;
    antiparticle_ = codenumber < 0;
    //if (! test_code()) {
    //  printf("%d %x\n", pdgcode_, pdgcode_);
    //  throw InvalidPDGCode("Invalid digits in PDG Code, they need to be Hex-encoded!");
    //}
  }

  /** Checks the integer for hex digits that are > 9.
   *
   * If one of the hex digits is not also a valid decimal digit,
   * something went wrong - maybe some user of this class forgot to
   * prefix the input with '0x' and thus passed 221 instead of 0x221.
   *
   * This routine only returns a boolean; it should be in the hands of
   * the calling function to raise this to an exception.
   *
   **/
  //inline bool test_code() {
  //  constexpr int bitmask = 0xf;
  //  // go through the hex digits:
  //  for (unsigned int bit = 0; bit < sizeof(int)*8; bit += 4) {
  //    // get the digit at hand: First, shift all bits to the small end,
  //    // so that the bits of interest are in the 4 least registers, then
  //    // use the bit mask to cut away the rest.
  //    int this_digit = (code_ >> bit) & bitmask;
  //    if (this_digit > 9 && bit < sizeof(int)*8 - 4) {
  //      return false;
  //    }
  //  }
  //  return true;
  //}

  inline int code() const {
    return (antiparticle_ << 31)
         | (n_    << 24)
         | (n_R_  << 20)
         | (n_L_  << 16)
         | (n_q1_ << 12)
         | (n_q2_ <<  8)
         | (n_q3_ <<  4)
         | (n_J_); 
  }

  /** The PDG code we're dealing with here.
   *
   * Everything else in this class is simpy smart accessors to this
   * integer.
   *
   * This integer is interpreted with hexadecimal digits, i.e., '545' is
   * interpreted as '0x221', which is what an eta-meson has (that 545 is
   * also its approximate mass in MeV is, I promise, a complete
   * coincidence).
   **/
   bool antiparticle_  : 1,
   // we don't use these bits.
                       : 3;
   std::uint32_t n_    : 4;
   std::uint32_t n_R_  : 4;
   std::uint32_t n_L_  : 4;
   std::uint32_t n_q1_ : 4;
   std::uint32_t n_q2_ : 4;
   std::uint32_t n_q3_ : 4;
   std::uint32_t n_J_  : 4;
 private:
};

} // namespace SMASH

#endif  // SRC_INCLUDE_PDGCODE_H_
