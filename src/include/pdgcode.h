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
#include <iostream>
#include <stdexcept>
#include <string>

namespace Smash {

/** PdgCode stores a Particle Data Group Particle Numbering Scheme
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
 * This class is simpy a collection of smart accessors to a 32 bit long
 * bitfield.
 *
 * The content is stored in hexadecimal digits, i.e., '545' is
 * interpreted as '0x221', which is what an eta-meson has (that 545 is
 * also its approximate mass in MeV is, I promise, a complete
 * coincidence).
 **/

class PdgCode {
 public:
  /// thrown for invalid inputs
  struct InvalidPdgCode : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };

  /// Standard initializer
  PdgCode() {
    try {
      set_fields(0xffffffffu);
    } catch (PdgCode::InvalidPdgCode) {
    }
  }
  /** Initialize using a string
   *
   * The string is interpreted as a hexadecimal number, i.e., @211@ is
   * interpreted as @0x211 = 529_{10}@.
   */
  PdgCode(const std::string& codestring) {
    set_from_string(codestring);
  }
  /** receive a signed integer and process it into a PDG Code. The sign
   * is taken as antiparticle boolean, while the absolute value of the
   * integer is used as hexdigits.
   */
  PdgCode(const int codenumber) {
    antiparticle_ = codenumber < 0;
    set_fields(static_cast<unsigned int>(std::abs(codenumber)));
  }
  /** receive an unsigned integer and process it into a PDG Code. The
   *  first bit is taken and used as antiparticle boolean.
   */
  PdgCode(const unsigned int abscode) {
    antiparticle_ = (abscode >> 31);
    set_fields(abscode);
  }

  /** Checks the integer for hex digits that are > 9.
   *
   * If one of the hex digits is not also a valid decimal digit,
   * something went wrong - maybe some user of this class forgot to
   * prefix the input with '0x' and thus passed 221 instead of 0x221.
   *
   * \return a bitmask indicating the offending digits. In the above
   * example, 221 = 0xd3, the second-to-last-digit is the offending one,
   * to the return value is 0b10 = 0x2 = 2.
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

  /** Dumps the bitfield into an unsigned integer. */
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

  /** Returns a signed integer with the PDG code in hexadecimal. */
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

  /// returns true if this is a baryon, antibaryon or meson.
  inline bool is_hadron() const {
    return (n_q3_ != 0 && n_q2_ != 0);
  }
  /// returns the baryon number of the particle.
  inline int baryon_number() const {
    if (! is_hadron() || n_q1_ == 0)
      return  0;
    return antiparticle_sign();
  }
  /** returns twice the isospin-3 component.
   *
   * This is calculated from the sum of heavy_quarkness of up and
   * down.
   */
  inline int isospin3() const {
    // heavy_quarkness(2) is the number of u quarks,
    // heavy_quarkness(1) is minus the number of d quarks.
    return heavy_quarkness(2)+heavy_quarkness(1);
  }
// TODO(baeuchle) This code needs special cases for things like η, p vs. Δ, Λ (vs. Σ)
//  /** returns twice the isospin vector length.
//   *
//   * This returns e.g. 2 for all pions and 3 for all Deltas. It is
//   * always positive.
//   */
//  inline unsigned int isospin_total() const {
//    // heavy_quarkness(2) is the number of u quarks,
//    // heavy_quarkness(1) is minus the number of d quarks.
//    return std::abs(heavy_quarkness(2))+std::abs(heavy_quarkness(1));
//  }
  /// returns the net number of \f$\bar s\f$ quarks.
  inline int strangeness() const {
    return heavy_quarkness(3);
  }
  /// returns the net number of \f$c\f$ quarks
  inline int charmness() const {
    return heavy_quarkness(4);
  }
  /// returns the net number of \f$\bar b\f$ quarks
  inline int bottomness() const {
    return heavy_quarkness(5);
  }
  /** returns the net number of (anti)quarks with given number
   *
   * Despite the name, it is also applicable to light quarks.
   *
   * \param quark PDG Code of quark: (1..8) = (d,u,s,c,b,t,b',t')
   * \return for u,c,t,t' quarks, it returns the net number of quarks (#quarks -
   * #antiquarks), while for d,s,b,b', it returns the net number of
   * antiquarks (#antiquarks - #quarks).
   *
   **/
  int heavy_quarkness(const int quark) const {
    // input sanitization: Only quark numbers 1 through 8 are allowed.
    if (quark < 1 || quark > 8)
      return 0;
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
  /** Returns the charge of the particle.
   *
   * The charge is calculated from the quark content (for hadrons) or
   * basically tabellized; currently leptons, neutrinos and the standard
   * model gauge bosons are known; unknown particles return a charge of
   * 0.
   **/
  int charge() const {
    if (is_hadron()) {
      // Q will accumulate 3*charge (please excuse the upper case. I
      // want to distinguish this from q which might be interpreted as
      // shorthand for "quark".)
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
  /** Returns twice the spin of a particle **/
  inline unsigned int spin() const {
    return n_J_-1;
  }
  /** Returns the spin degeneracy \f$2s + 1\f$ of a particle **/
  inline unsigned int spin_degeneracy() const {
    return n_J_;
  }
  /// returns -1 for antiparticles and +1 for particles.
  inline int antiparticle_sign() const {
    return (antiparticle_ ? -1 : +1);
  }

  /// istream >> PdgCode assigns the PDG Code from an istream.
  friend std::istream& operator>>(std::istream& is, PdgCode& code);

 private:
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

  /** extract digits from a character. */
  inline std::uint32_t get_digit_from_char(const char inp) const {
    // atoi's behaviour for invalid input is undefined. I don't like
    // that.

    // this checks if the first four digits are 0011 (as they should be
    // for ASCII digits).
    if ((inp & 0xf0) ^ 0x30) {
      throw InvalidPdgCode("Invalid character " + std::to_string(inp)
                         + " found.\n");
    }
    // the last four digits are the number; they should not be > 9
    // (i.e., one of [:;<=>?])
    if ((inp & 0x0f) > 9) {
      throw InvalidPdgCode("Invalid digit " + std::to_string(inp)
                         + " found.\n");
    }
    // now that we've checked that the first bits are correct and the
    // last bits are a number, we can return the last bits.
    return (inp & 0x0f);
  }

  // takes a string and sets the fields.
  inline void set_from_string(const std::string& codestring) {
    antiparticle_ = false;
    n_ = n_R_ = n_L_ = n_q1_ = n_q2_ = n_q3_ = n_J_ = 0;
    size_t length = codestring.size();
    if (length < 1) {
      throw InvalidPdgCode("Empty string does not contain PDG Code\n");
    }
    int c = 0;
    // look at current character; if it is a + or minus sign, read it
    // and advance to next char.
    if (codestring[c] == '-') {
      antiparticle_ = true;
      c++;
    } else if (codestring[c] == '+') {
      antiparticle_ = false;
      c++;
    }
    // save if the first character was a sign:
    unsigned int sign = c;
    // codestring shouldn't be longer than 7 + sign.
    if (length > 7+sign) {
      throw InvalidPdgCode("String \"" + codestring +
                           "\" too long for PDG Code\n");
    }
    // codestring has 7 digits? 7th from last goes in n_.
    if (length > 6+sign) {
      n_ = get_digit_from_char(codestring[c++]);
    }
    // it has 6 or 7 digits? 6th from last is n_R_.
    if (length > 5+sign) {
      n_R_ = get_digit_from_char(codestring[c++]);
    }
    // 5th from last is n_L_.
    if (length > 4+sign) {
      n_L_ = get_digit_from_char(codestring[c++]);
    }
    // 4th from last is n_q1_.
    if (length > 3+sign) {
      n_q1_ = get_digit_from_char(codestring[c++]);
    }
    // 3rd from last is n_q2_.
    if (length > 2+sign) {
      n_q2_ = get_digit_from_char(codestring[c++]);
    }
    // next to last is n_q3_.
    if (length > 1+sign) {
      n_q3_ = get_digit_from_char(codestring[c++]);
    }
    // last digit is the spin degeneracy.
    if (length > sign) {
      n_J_ = get_digit_from_char(codestring[c++]);
    } else {
      throw InvalidPdgCode("String \"" + codestring +
                 "\" only consists of a sign, that is no valid PDG Code\n");
    }
  }

  /** Sets the bitfield from an unsigned integer. Usually called from
   * the constructors.
   **/
  inline void set_fields(unsigned int abscode) {
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
    int test = test_code();
    if (test > 0) {
      throw InvalidPdgCode("Invalid digits " + std::to_string(test) + 
                           " in PDG Code " + std::to_string(code()));
    }
  }
};

std::istream& operator>>(std::istream& is, PdgCode& code) {
  std::string codestring("");
  // discard any whitespace at beginning:
  while (is.peek() == ' ') {
    is.get();
  }
  // read sign if there is one:
  if (is.peek() == '+' || is.peek() == '-') {
    codestring += is.get();
  }
  // read a maximum of 7 characters
  for (int c = 0; c < 7; c++) {
    // look into the string. Is it a valid character?
    try {
      code.get_digit_from_char(is.peek());
    } catch (PdgCode::InvalidPdgCode) {
      // if not, end the loop.
      break;
    }
    // read one character from is:
    // char * s = new char[1];
    // is.read(s, 1);
    codestring += is.get();
  }
  try {
    // set the fields from the string:
    code.set_from_string(codestring);
  } catch (PdgCode::InvalidPdgCode) {
    is.setstate(std::ios::failbit);
  }
  return is;
}

} // namespace SMASH

#endif  // SRC_INCLUDE_PDGCODE_H_
