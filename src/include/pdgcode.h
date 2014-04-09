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
    unsigned int abscode = std::abs(codenumber);
    antiparticle_ = codenumber < 0;
    set_fields(abscode);
    int test = test_code();
    if (test > 0) {
      throw InvalidPDGCode("Invalid digits in PDG Code " + std::to_string(code()));
     // throw InvalidPDGCode(std::string("Invalid digits in PDG Code ")
     //                    + code() + std::string(" they need to be Hex-encoded!"));
    }
  }
  PDGCode(unsigned int abscode) {
    set_fields(abscode);
    abscode >>= 31;
    antiparticle_ = abscode;
    int test = test_code();
    if (test > 0) {
      throw InvalidPDGCode("Invalid digits in PDG Code " + std::to_string(code()));
     // throw InvalidPDGCode(std::string("Invalid digits in PDG Code ")
     //                    + code() + std::string(" they need to be Hex-encoded!"));
    }
  }

  inline void set_fields(unsigned int& abscode) {
    constexpr unsigned int bitmask = 0xf;
    n_J_ = abscode & bitmask;
    abscode >>= 4;
    n_q3_ = abscode & bitmask;
    abscode >>= 4;
    n_q2_ = abscode & bitmask;
    abscode >>= 4;
    n_q1_ = abscode & bitmask;
    abscode >>= 4;
    n_L_  = abscode & bitmask;
    abscode >>= 4;
    n_R_  = abscode & bitmask;
    abscode >>= 4;
    n_    = abscode & bitmask;
  }

  /** Checks the integer for hex digits that are > 9.
   *
   * If one of the hex digits is not also a valid decimal digit,
   * something went wrong - maybe some user of this class forgot to
   * prefix the input with '0x' and thus passed 221 instead of 0x221.
   *
   * \return fail a bitmask indicating the offending digits. In the
   * above example, 221 = 0xd3, the second-to-last-digit is the
   * offending one, to the return value is 0b10 = 0x2 = 2.
   *
   **/
  inline int test_code() {
    int fail = 0;
    if (n_    > 9) fail |= 1<<6;
    if (n_R_  > 9) fail |= 1<<5;
    if (n_L_  > 9) fail |= 1<<4;
    if (n_q1_ > 9) fail |= 1<<3;
    if (n_q2_ > 9) fail |= 1<<2;
    if (n_q3_ > 9) fail |= 1<<1;
    if (n_J_  > 9) fail |= 1;
    return fail;
  }

  inline unsigned int dump() const {
    return (antiparticle_ << 31)
         | (n_    << 24)
         | (n_R_  << 20)
         | (n_L_  << 16)
         | (n_q1_ << 12)
         | (n_q2_ <<  8)
         | (n_q3_ <<  4)
         | (n_J_); 
  }

  inline int code() const {
    return (antiparticle_ ? -1 : 1)
         * ( (n_    << 24)
           | (n_R_  << 20)
           | (n_L_  << 16)
           | (n_q1_ << 12)
           | (n_q2_ <<  8)
           | (n_q3_ <<  4)
           | (n_J_)); 
  }

  /// returns true if this is a baryon OR antibaryon!.
  /// returns true if \f$n_{q_1} \ne 0\f$ and \f$n_{q_2} \ne 0\f$.
  inline bool is_hadron() const {
    return (n_q3_ != 0 && n_q2_ != 0);
  }
  /// returns true if this is a meson.
  inline int baryon_number() const {
    if (! is_hadron() || n_q1_ == 0)
      return  0;
    return antiparticle_sign();
  }
  /// returns twice the isospin-3 component.
  inline int isospin3() const {
    // heavy_quarkness(2) is the number of u quarks,
    // heavy_quarkness(1) is minus the number of d quarks.
    return heavy_quarkness(2)+heavy_quarkness(1);
  }
  inline int strangeness() const {
    return heavy_quarkness(3);
  }
  inline int charmness() const {
    return heavy_quarkness(4);
  }
  inline int bottomness() const {
    return heavy_quarkness(5);
  }
  /** this looks for the quantum number associated with quark number
   * quark.
   *
   * Despite the name, it is also applicable to light quarks.
   *
   **/
  int heavy_quarkness(const int quark) const {
    // non-hadrons and those that have none of this quark type: 0.
    if (! is_hadron() || (n_q1_ != quark && n_q2_ != quark && n_q3_ != quark))
      return 0;
    // baryons: count quarks.
    if (baryon_number() != 0) {
      // u,s,b quarks get negative *ness; d,c,t have positive *ness:
      int sign = (quark%2 == 0) ? +1 : -1;
      // and for anti-baryons, the sign changes:
      sign *= antiparticle_sign();
      return sign*((n_q1_ == quark) + (n_q2_ == quark) + (n_q3_ == quark));
    }
    // mesons.
    // quarkonium state? Not open heavy-quarkness.
    if (n_q3_ == quark && n_q2_ == quark)
      return 0;
    // this has covered all the easy stuff
    // get the "other" quark. (We know this must exist, since they are
    // not both the right one and one of them is the right one).
    int otherquark = (n_q2_ == quark) ? n_q3_ : n_q2_;
    // "our" quark is the heavier one: 1 for particles
    if (otherquark < quark) {
      return antiparticle_sign();
    }
    // ours is the lighter: If the heavier quark has the same SIGN, we
    // get a -1, else +1.
    if (otherquark%2 == quark%2)
      return -1*antiparticle_sign();
    return antiparticle_sign();
  }
  int charge() const {
    if (is_hadron()) {
      // this will accumulate 3*charge
      int Q = 0;
      // this loops over d,u,s,c,b,t,b',t' quarks (the latter three can
      // be safely ignored, but I don't think this will be a bottle
      // neck.
      for (int i = 1; i < 9; i++) {
        // the appropriate sign is already in heavy_quarkness; now
        // u,c,t,t' quarks have charge = 2/3 e, while d,s,b,b' quarks
        // have -1/3 e.
        Q += (i%2 == 0 ? 2 : 1)*heavy_quarkness(i);
      }
      return Q/3;
    }
    // non-hadron:
    // Leptons: 11, 13, 15, 17 are e, μ, τ, τ' and have a charge -1,
    // while    12, 14, 16, 18 are the neutrinos that have no charge.
    if (n_q3_ == 1)
      return -1*(n_J_%2)*antiparticle_sign();
    // Bosons: 24 is the W+, all else is uncharged.
    if (n_q3_ == 2 && n_J_ == 4)
      return antiparticle_sign();
    // default (this includes all other Bosons) is 0.
    return 0;
  }
  inline int antiparticle_sign() const {
    return (antiparticle_ ? -1 : +1);
  }

 private:
  /** The PDG code we're dealing with here.
   *
   * Everything else in this class is simpy smart accessors to this
   * bitfield.
   *
   * The number is interpreted with hexadecimal digits, i.e., '545' is
   * interpreted as '0x221', which is what an eta-meson has (that 545 is
   * also its approximate mass in MeV is, I promise, a complete
   * coincidence).
   *
   * These digits are stored in the bit fields.
   **/
   /// first bit: stores the sign.
   bool antiparticle_  : 1,
   // we don't use these bits.
                       : 3;
   /// first field: "counter"
   std::uint32_t n_    : 4;
   /// "radial excitation"
   std::uint32_t n_R_  : 4;
   /// "angular momentum"
   std::uint32_t n_L_  : 4;
   /// first quark field. 0 for mesons.
   std::uint32_t n_q1_ : 4;
   /// second quark field
   std::uint32_t n_q2_ : 4;
   /// third quark field
   std::uint32_t n_q3_ : 4;
   /// spin quantum number \f$n_J = 2 J + 1\f$.
   std::uint32_t n_J_  : 4;
};

} // namespace SMASH

#endif  // SRC_INCLUDE_PDGCODE_H_
