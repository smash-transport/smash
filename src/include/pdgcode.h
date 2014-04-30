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
 * \see http://pdg.lbl.gov/2013/reviews/rpp2012-rev-monte-carlo-numbering.pdf
 *
 * Usage:
 * ------
 * \code
 * #include "include/pdgcode.h"
 *
 * // initalize with an integer: make sure it is hex-encoded!
 * PdgCode pi_plus(0x211);
 * // you can also initalize from a string:
 * PdgCode pi_minus("-211");
 * // initialize a PDG Code that knows it is not set yet:
 * PdgCode other_particle();
 * // this is true:
 * if (other_particle == PdgCode::invalid()) {
 *   printf("Invalid particle! Please enter PDG Code: ");
 *   // fill from stringstream:
 *   std::cin >> other_particle;
 * }
 * // check if particle is in the same multiplet as the pions:
 * if (pi_plus.multiplet() == other_particle.multiplet()) {
 *   printf("The other particle is a 0^-- meson.\n");
 * }
 * // is this a Kaon?
 * if (other_particle.code() == 0x311) {
 *   printf("The particle is a K plus\n");
 * }
 * // what baryon number does the particle have?
 * printf("The particle has a baryon number of %d\n",
 *                                           other_particle.baryon_number());
 * \endcode
 *
 * This class contains a collection of smart accessors to the PDG code
 * so that quantum numbers etc can easily be read off.
 *
 * Internals
 * ---------
 *
 * The content is stored in hexadecimal digits, i.e., the number '545'
 * is interpreted as '0x221', i.e., an eta-meson. To check if a given
 * particle is of a given type, make sure that you give the type in hex
 * digits as well (see example above).
 *
 * The reason for that is that the concept of PdgCodes, especially for
 * Hadrons, is not one of wholesale numbers, but one of concatenated
 * digits. Using hexadecimally interpreted digits makes it numerically
 * very easy to access the separate digits (there's no arithmetic
 * involved with successive divisions by 10 and taking the remainder
 * etc.).
 *
 * Limitations:
 * ------------
 *
 * The code is tuned to non-colored objects at the moment. That means
 * that colored objects (Diquarks and Quarks) are not easily useable
 * with this class; the behaviour of functions baryon_number, charge,
 * is_hadron etc. is undefined. (This is mostly because these things are
 * not well-defined, and/or because the charge and baryon number is not
 * an integer anymore.)
 *
 * Also, tetra- and pentaquarks as well as compound objects (small or
 * large nuclei) cannot be represented; that, though, is a problem of
 * the PDG Numbering Scheme rather than of this class.
 *
 **/

class PdgCode {
 public:
  /// thrown for invalid inputs
  struct InvalidPdgCode : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };

  /****************************************************************************
   *                                                                          *
   * First, the constructors                                                  *
   *                                                                          *
   ****************************************************************************/
  /// Standard initializer
  PdgCode() : dump_(0x0) {}
  /** Initialize using a string
   *
   * The string is interpreted as a hexadecimal number, i.e., \c 211 is
   * interpreted as \c 0x211 = \f$529_{10}\f$.
   */
  PdgCode(const std::string& codestring) {
    set_from_string(codestring);
  }
  /** receive a signed integer and process it into a PDG Code. The sign
   * is taken as antiparticle boolean, while the absolute value of the
   * integer is used as hexdigits.
   *
   * \param codenumber The number 0x221 is interpreted as an η meson,
   * -0x211 is a "charged pi antiparticle", i.e., a \f$\pi^-\f$.
   */
  PdgCode(std::int32_t codenumber) : dump_(0x0)  {
    digits_.antiparticle_ = false;
    if (codenumber < 0) {
      digits_.antiparticle_ = true;
      codenumber = -codenumber;
    }
    set_fields(codenumber);
  }
  /** receive an unsigned integer and process it into a PDG Code. The
   *  first bit is taken and used as antiparticle boolean.
   */
  PdgCode(const std::uint32_t abscode) : dump_(0x0) {
    // use the first bit for the antiparticle_ boolean.
    digits_.antiparticle_ = ((abscode & 0x80000000u) != 0);
    set_fields(abscode);
  }

  /****************************************************************************
   *                                                                          *
   * test function and export functions                                       *
   *                                                                          *
   ****************************************************************************/
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
    if (digits_.n_    > 9) { fail |= 1<<6; }
    if (digits_.n_R_  > 9) { fail |= 1<<5; }
    if (digits_.n_L_  > 9) { fail |= 1<<4; }
    if (digits_.n_q1_ > 9) { fail |= 1<<3; }
    if (digits_.n_q2_ > 9) { fail |= 1<<2; }
    if (digits_.n_q3_ > 9) { fail |= 1<<1; }
    if (digits_.n_J_  > 9) { fail |= 1; }
    return fail;
  }

  /** Dumps the bitfield into an unsigned integer. */
  inline std::uint32_t dump() const {
    // this cuts the three unused bits.
    return (dump_ & 0x8fffffff);
  }

  /** Returns a signed integer with the PDG code in hexadecimal. */
  inline std::int32_t code() const {
    return antiparticle_sign() * (dump_ & 0x0fffffff);
  }

  /// returns a C++ string from the PDG Code.
  inline std::string string() const {
    char hexstring[8];
    if (digits_.antiparticle_) {
      snprintf(hexstring, 8, "-%x", std::abs(code()));
    } else {
      snprintf(hexstring, 8, "%x", code());
    }
    return std::string(hexstring);
  }

  /****************************************************************************
   *                                                                          *
   * accessors of various properties                                          *
   *                                                                          *
   ****************************************************************************/
  /// returns true if this is a baryon, antibaryon or meson.
  inline bool is_hadron() const {
    return (digits_.n_q3_ != 0 && digits_.n_q2_ != 0);
  }
  /// returns the baryon number of the particle.
  inline int baryon_number() const {
    if (! is_hadron() || digits_.n_q1_ == 0) {
      return  0;
    }
    return antiparticle_sign();
  }
  /** returns twice the isospin-3 component \f$I_3\f$.
   *
   * This is calculated from the sum of net_quark_number of up and down.
   */
  inline int isospin3() const {
    // net_quark_number(2) is the number of u quarks,
    // net_quark_number(1) is the number of d quarks.
    return net_quark_number(2)-net_quark_number(1);
  }
  /** returns twice the isospin vector length \f$I\f$.
   *
   * This returns e.g. 2 for all pions and 3 for all Deltas. It is
   * always positive.
   */
  unsigned int isospin_total() const;
  /** returns the net number of \f$\bar s\f$ quarks.
   *
   * For particles with one strange quark, -1 is returned.
   */
  inline int strangeness() const {
    return -net_quark_number(3);
  }
  /** returns the net number of \f$c\f$ quarks
   *
   * For particles with one charm quark, +1 is returned.
   */
  inline int charmness() const {
    return +net_quark_number(4);
  }
  /** returns the net number of \f$\bar b\f$ quarks
   *
   * For particles with one bottom quark, -1 is returned.
   */
  inline int bottomness() const {
    return -net_quark_number(5);
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
        // u,c,t,t' quarks have charge = 2/3 e, while d,s,b,b' quarks
        // have -1/3 e. The antiparticle sign s already in net_quark_number.
        Q += (i % 2 == 0 ? 2 : -1) * net_quark_number(i);
      }
      return Q / 3;
    }
    // non-hadron:
    // Leptons: 11, 13, 15, 17 are e, μ, τ, τ' and have a charge -1,
    // while    12, 14, 16, 18 are the neutrinos that have no charge.
    if (digits_.n_q3_ == 1) {
      return -1 * (digits_.n_J_ % 2) * antiparticle_sign();
    }
    // Bosons: 24 is the W+, all else is uncharged.
    // we ignore the first digits so that this also finds strange gauge
    // boson "resonances" (in particular, \f$\tilde \chi_1^+\f$ with PDG
    // Code 1000024).
    if ((dump_ & 0x0000ffff) == 0x24) {
      return antiparticle_sign();
    }
    // default (this includes all other Bosons) is 0.
    return 0;
  }
  /** Returns twice the spin of a particle.
   *
   * The code is good for hadrons, leptons and spin-1-bosons. It returns
   * 2 (meaning spin=1) for the Higgs, though.
   */
  inline unsigned int spin() const {
    if (is_hadron()) {
      return digits_.n_J_ - 1;
    }
    // this assumes that we only have white particles (no single
    // quarks): Electroweak fermions have 11-17, so the
    // second-to-last-digit is the spin. The same for the Bosons: they
    // have 21-29 and 2spin = 2 (this fails for the Higgs).
    return digits_.n_q3_;
  }
  /** Returns the spin degeneracy \f$2s + 1\f$ of a particle **/
  inline unsigned int spin_degeneracy() const {
    if (is_hadron()) {
      return digits_.n_J_;
    }
    return spin() + 1;
  }
  /// returns -1 for antiparticles and +1 for particles.
  inline int antiparticle_sign() const {
    return (digits_.antiparticle_ ? -1 : +1);
  }
  /// returns an integer with only the quark numbers set.
  inline std::int32_t quarks() const {
    if (! is_hadron()) {
      return 0;
    }
    return chunks_.quarks_;
  }
  /** Returns an identifier for the SU(N) multiplet of this PDG Code.
   *
   * This can be used to compare if two particles are in the same SU(N)
   * multiplet, i.e., p, Λ, Σ and \f$\Xi_{cb}\f$ are in the same
   * multiplet, as are Δ, Ω, \f$\Omega_{ccc}\f$ and ρ, ω,
   * \f$D^\ast_s\f$.
   *
   * The antiparticle sign is ignored for mesons: Both \f$K^+\f$ and
   * \f$K^-\f$ are in the same multiplet; the same applies for charmed
   * and anti-charmed mesons.
   *
   * The exact format of this is subject to change; only use this to
   * check if two particles are in the same multiplet!
   */
  inline std::int32_t multiplet() const {
    if (! is_hadron()) {
      return 0;
    }
    // the multiplet code is the antiparticle_*( baryon number +
    // excitation_ + n_J_) [the "+" being a concatenation here].  Baryon
    // number and sign are added later.
    std::int32_t multiplet_code = ((chunks_.excitation_ << 4) | digits_.n_J_);
    // if the particle is in a meson multiplet, there are no signs.
    if (baryon_number() == 0) {
      return multiplet_code;
    }
    // else, return the sign and the baryon number, too.
    return antiparticle_sign() *
                  (multiplet_code | ((baryon_number() & 1) << 16));
  }
  /** Returns an identifier for the Isospin-multiplet of this PDG Code.
   *
   * This can be used to compare if two particles are in the same
   * isospin multiplet.  For non-hadrons, this returns 0, and for
   * hadrons, it returns the PDG Code of the similar particle that has
   * all up quarks replaced with down quarks.
   *
   * To distinguish η mesons from \f$\pi^0\f$ and Λ from \f$\Sigma^0\f$, Isopin
   * = 0 states return the complete code.
   *
   * The antiparticle sign is ignored for pi-like mesons, so that \f$\pi^+\f$
   * and \f$\pi^-\f$ are in the same multiplet. The sign is used for all other
   * mesons and baryons, so that \f$K^-\f$ and \f$K^+\f$ are not in the same
   * multiplet, and neither \f$p\f$ and \f$\bar p\f$.
   *
   */
  inline std::int32_t iso_multiplet() const {
    if (! is_hadron()) {
      return 0;
    }
    // the η meson is NOT in the same multiplet as the π, and Λ and Σ.
    // We return the complete code for the I=0 state.
    if (isospin_total() == 0) {
      return code();
    }
    // this is the code with all u quarks changed to d quarks. Note that it
    // doesn't matter that we change e.g. a proton to an (ddd)-state, because
    // we can distinguish nucleons and Δ resonances by the other numbers in the
    // new scheme (it's gotta be good for something!)
    std::int32_t multiplet_code = ( (chunks_.excitation_  << 16)
                          | ((digits_.n_q1_ == 2 ? 1 : digits_.n_q1_) << 12)
                          | ((digits_.n_q2_ == 2 ? 1 : digits_.n_q2_) <<  8)
                          | ((digits_.n_q3_ == 2 ? 1 : digits_.n_q3_) <<  4)
                          | (digits_.n_J_));
    // if we have pion-like particles, return the above code (discard
    // antiparticle_sign)
    if ((multiplet_code & 0x0000fff0) == 0x110) {
      return multiplet_code;
    }
    // else, the sign is important!
    return antiparticle_sign()*multiplet_code;
  }

  /****************************************************************************
   *                                                                          *
   * operations with more than one PDG Code                                   *
   *                                                                          *
   ****************************************************************************/
  /** sorts PDG Codes according to their numeric value.
   *
   * This is used by std::map
   */
  inline bool operator<(const PdgCode rhs) const {
    return (code() < rhs.code());
  }
  /// returns if the codes are equal
  inline bool operator==(const PdgCode rhs) const {
    return (code() == rhs.code());
  }
  /// returns if the codes are not equal.
  inline bool operator!=(const PdgCode rhs) const {
    return ! (*this == rhs);
  }
  /// returns if the code of rhs is the inverse of this one.
  inline bool is_antiparticle_of(const PdgCode rhs) const {
    return code() == -rhs.code();
  }

  /// istream >> PdgCode assigns the PDG Code from an istream.
  friend std::istream& operator>>(std::istream& is, PdgCode& code);

  /** PdgCode 0x0 is guaranteed not to be valid by the PDG standard, but
   * it passes all tests here, so we can use it to show some code is not
   * yet set.
   */
  static PdgCode invalid() { return PdgCode(0x0); }

 private:
#if !defined(__GNUC__) || !defined(__x86_64__)
#error "Please determine the correct bit-field order for your target/compiler"
#endif
  /** the union holds the data; either as a single integer dump_, as a
   * single-digit bitfield digits_ or as a multiple-digits bitfield
   * chunks_.
   */
  union {
    /** the single digits collection of the code. Here, every PDG code
     * digits is directly accessible. */
    struct {
#if defined(__GNUC__) || defined(__x86_64__) || defined(DOXYGEN)
      /// spin quantum number \f$n_J = 2 J + 1\f$.
      std::uint32_t n_J_  : 4;
      /// third quark field
      std::uint32_t n_q3_ : 4;
      /// second quark field
      std::uint32_t n_q2_ : 4;
      /// first quark field. 0 for mesons.
      std::uint32_t n_q1_ : 4;
      /// "angular momentum"
      std::uint32_t n_L_  : 4;
      /// "radial excitation"
      std::uint32_t n_R_  : 4;
      /// first field: "counter"
      std::uint32_t n_    : 4, :3;
      /// first bit: stores the sign.
      bool antiparticle_  : 1;
#else  // reverse ordering
      bool antiparticle_  : 1, :3;
      std::uint32_t n_    : 4;
      std::uint32_t n_R_  : 4;
      std::uint32_t n_L_  : 4;
      std::uint32_t n_q1_ : 4;
      std::uint32_t n_q2_ : 4;
      std::uint32_t n_q3_ : 4;
      std::uint32_t n_J_  : 4;
#endif  // __GNUC__ or __x86_64__
    } digits_;
    /** the bitfield dumped into a single integer. Please note that the
     * 2nd, 3rd and 4th highest bits are possibly undefined.
     **/
    std::uint32_t dump_;
    /** chunk collection: here, the chunks with \f$nn_Rn_L\f$ and
     * \f$n_{q_1}n_{q_2}n_{q_3}\f$ are directly accessible.
     */
    struct {
#if defined(__GNUC__) || defined(__x86_64__) || defined(DOXYGEN)
      std::uint32_t             :  4;
      /// the quark digits n_q{1,2,3}_
      std::uint32_t quarks_     : 12;
      /// the excitation digits n_, n_R_, n_L_
      std::uint32_t excitation_ : 12, : 4;
#else
      std::uint32_t : 4, excitation_ : 12;
      std::uint32_t quarks_     : 12, : 4;
#endif  // __GNUC__ or __x86_64__
    } chunks_;
  };

  /** returns the net number of quarks with given flavour number
   *
   * \param quark PDG Code of quark: (1..8) = (d,u,s,c,b,t,b',t')
   * \return for the net number of quarks (\#quarks - \#antiquarks)
   *
   * For public use, see strangeness(), charmness(), bottomness() and
   * isospin3().
   **/
  int net_quark_number(const int quark) const;
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

  /// takes a string and sets the fields.
  inline void set_from_string(const std::string& codestring) {
    digits_.antiparticle_ = false;
    digits_.n_ = digits_.n_R_ = digits_.n_L_ = digits_.n_q1_ = digits_.n_q2_ =
                                digits_.n_q3_ = digits_.n_J_ = 0;
    size_t length = codestring.size();
    if (length < 1) {
      throw InvalidPdgCode("Empty string does not contain PDG Code\n");
    }
    int c = 0;
    // look at current character; if it is a + or minus sign, read it
    // and advance to next char.
    if (codestring[c] == '-') {
      digits_.antiparticle_ = true;
      ++c;
    } else if (codestring[c] == '+') {
      digits_.antiparticle_ = false;
      ++c;
    }
    // save if the first character was a sign:
    unsigned int sign = c;
    // codestring shouldn't be longer than 7 + sign.
    if (length > 7 + sign) {
      throw InvalidPdgCode("String \"" + codestring +
                           "\" too long for PDG Code\n");
    }
    // please note that in what follows, we actually need c++, not ++c.
    // codestring has 7 digits? 7th from last goes in n_.
    if (length > 6 + sign) {
      digits_.n_ = get_digit_from_char(codestring[c++]);
    }
    // it has 6 or 7 digits? 6th from last is n_R_.
    if (length > 5 + sign) {
      digits_.n_R_ = get_digit_from_char(codestring[c++]);
    }
    // 5th from last is n_L_.
    if (length > 4 + sign) {
      digits_.n_L_ = get_digit_from_char(codestring[c++]);
    }
    // 4th from last is n_q1_.
    if (length > 3 + sign) {
      digits_.n_q1_ = get_digit_from_char(codestring[c++]);
    }
    // 3rd from last is n_q2_.
    if (length > 2 + sign) {
      digits_.n_q2_ = get_digit_from_char(codestring[c++]);
    }
    // next to last is n_q3_.
    if (length > 1 + sign) {
      digits_.n_q3_ = get_digit_from_char(codestring[c++]);
    }
    // last digit is the spin degeneracy.
    if (length > sign) {
      digits_.n_J_ = get_digit_from_char(codestring[c++]);
    } else {
      throw InvalidPdgCode("String \"" + codestring +
                 "\" only consists of a sign, that is no valid PDG Code\n");
    }
  }

  /** Sets the bitfield from an unsigned integer. Usually called from
   * the constructors.
   **/
  inline void set_fields(std::uint32_t abscode) {
    // "dump_ =" overwrites antiparticle_, but this needs to have been set
    // already, so we carry it around the assignment.
    bool ap = digits_.antiparticle_;
    dump_ = abscode & 0x0fffffff;
    digits_.antiparticle_ = ap;
    int test = test_code();
    if (test > 0) {
      throw InvalidPdgCode("Invalid digits " + std::to_string(test) +
                           " in PDG Code " + std::to_string(code()));
    }
  }

};

std::istream& operator>>(std::istream& is, PdgCode& code);

}  // namespace SMASH

#endif  // SRC_INCLUDE_PDGCODE_H_
