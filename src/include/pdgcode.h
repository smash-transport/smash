/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_PDGCODE_H_
#define SRC_INCLUDE_PDGCODE_H_

#include <algorithm>
#include <array>
#include <cstdlib>
#include <iosfwd>
#include <sstream>
#include <stdexcept>
#include <string>

#include "pdgcode_constants.h"

namespace Smash {

/**
 * \ingroup data
 *
 * PdgCode stores a Particle Data Group Particle Numbering Scheme
 * particle type number.
 *
 * \see http://pdg.lbl.gov/2014/reviews/rpp2014-rev-monte-carlo-numbering.pdf
 *
 * Usage:
 * ------
 * \code
 * #include "include/pdgcode.h"
 *
 * // initialize with an integer: make sure it is hex-encoded!
 * PdgCode pi_plus(0x211);
 * // you can also initialize from a string:
 * PdgCode pi_minus("-211");
 * // initialize a PDG Code that knows it is not set yet:
 * PdgCode other_particle();
 * // this is true:
 * if (other_particle == PdgCode::invalid()) {
 *   printf("Invalid particle! Please enter PDG Code: ");
 *   // fill from stringstream:
 *   std::cin >> other_particle;
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
  /// \ingroup exception
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
  explicit PdgCode(const std::string& codestring) {
    set_from_string(codestring);
  }

  /** receive a signed integer and process it into a PDG Code. The sign
   * is taken as antiparticle boolean, while the absolute value of the
   * integer is used as hexdigits.
   *
   * \param codenumber The number 0x221 is interpreted as an η meson,
   * -0x211 is a "charged pi antiparticle", i.e., a \f$\pi^-\f$.
   */
  PdgCode(std::int32_t codenumber) : dump_(0x0) {  // NOLINT(runtime/explicit)
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
  explicit PdgCode(const std::uint32_t abscode) : dump_(0x0) {
    // use the first bit for the antiparticle_ boolean.
    digits_.antiparticle_ = ((abscode & 0x80000000u) != 0);
    set_fields(abscode);
  }

  /****************************************************************************
   *                                                                          *
   * test function and export functions                                       *
   *                                                                          *
   ****************************************************************************/
  /** Checks the integer for invalid hex digits.
   *
   * Usually all digits are at least <= 9. The n_q digits are even <= 6
   * (because there are only six quarks).
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
  inline int test_code() const {
    int fail = 0;
    if (digits_.n_ > 9) {
      fail |= 1 << 6;
    }
    if (digits_.n_R_ > 9) {
      fail |= 1 << 5;
    }
    if (digits_.n_L_ > 9) {
      fail |= 1 << 4;
    }
    if (digits_.n_q1_ > 6) {
      fail |= 1 << 3;
    }
    if (digits_.n_q2_ > 6) {
      fail |= 1 << 2;
    }
    if (digits_.n_q3_ > 6) {
      fail |= 1 << 1;
    }
    if (digits_.n_J_ > 9) {
      fail |= 1;
    }
    return fail;
  }

  /// Do all sorts of validity checks.
  void check() const {
    // n_J must be odd for mesons and even for baryons (and cannot be zero)
    if (is_hadron()) {
      if (baryon_number() == 0) {
        // mesons: special cases K0_L=0x130 and K0_S=0x310
        if ((digits_.n_J_ % 2 == 0) && dump() != 0x130 && dump() != 0x310) {
          throw InvalidPdgCode("Invalid PDG code " + string() +
                               " (meson with even n_J)");
        }
      } else {
        if ((digits_.n_J_ % 2 != 0) || digits_.n_J_ == 0) {
          throw InvalidPdgCode("Invalid PDG code " + string() +
                               " (baryon with odd n_J)");
        }
      }
    } else {
      if (digits_.n_J_ == 0 && dump() != 0x0) {
        throw InvalidPdgCode("Invalid PDG code " + string() + " (n_J==0)");
      }
    }
    /* The antiparticle flag only makes sense for particle types
     * that have an antiparticle. */
    if (digits_.antiparticle_ && !has_antiparticle()) {
      throw InvalidPdgCode("Invalid PDG code " + string() +
                           " (cannot be negative)");
    }
  }

  /** Dumps the bitfield into an unsigned integer. */
  inline std::uint32_t dump() const {
    // this cuts the three unused bits.
    return (dump_ & 0x8fffffff);
  }

  /** Returns a signed integer with the PDG code in hexadecimal. */
  inline std::int32_t code() const { return antiparticle_sign() * ucode(); }

  /// returns a C++ string from the PDG Code.
  inline std::string string() const {
    std::stringstream ss;
    if (digits_.antiparticle_) {
      ss << "-";
    }
    ss << std::hex << ucode();
    return ss.str();
  }

  /// Construct the antiparticle to a given PDG code.
  PdgCode get_antiparticle() const {
    // TODO(mkretz): more efficient implementation
    return PdgCode(-code());
  }

  /// Construct PDG code from decimal number
  static PdgCode from_decimal(const int pdgcode_decimal) {
    int a = pdgcode_decimal;
    int hex_pdg = 0, tmp = 1;
    while (a) {
      hex_pdg += (a % 10) * tmp;
      tmp *= 16;
      a = a / 10;
    }
    return PdgCode(hex_pdg);
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
  /// returns true if this is a lepton.
  inline bool is_lepton() const {
    return (digits_.n_q1_ == 0 && digits_.n_q2_ == 0 && digits_.n_q3_ == 1);
  }
  /// returns the baryon number of the particle.
  inline int baryon_number() const {
    if (!is_hadron() || digits_.n_q1_ == 0) {
      return 0;
    }
    return antiparticle_sign();
  }
  /// Returns whether this PDG code identifies a baryon.
  inline bool is_baryon() const { return is_hadron() && digits_.n_q1_ != 0; }

  /// Returns whether this PDG code identifies a meson.
  inline bool is_meson() const { return is_hadron() && digits_.n_q1_ == 0; }

  /// Is this a nucleon/anti-nucleon (p, n, -p, -n)?
  inline bool is_nucleon() const {
    const auto abs_code = std::abs(code());
    return (abs_code == pdg::p || abs_code == pdg::n);
  }

  /// Is this a N*(1535) (+/0)?
  inline bool is_Nstar1535() const {
    const auto abs_code = std::abs(code());
    return (abs_code == pdg::N1535_p || abs_code == pdg::N1535_z);
  }

  /// Is this a Delta(1232) (with anti-delta)?
  inline bool is_Delta() const {
    const auto abs_code = std::abs(code());
    return (abs_code == pdg::Delta_pp || abs_code == pdg::Delta_p ||
            abs_code == pdg::Delta_z || abs_code == pdg::Delta_m);
  }

  /// Is this a hyperon (Lambda, Sigma, Xi, Omega)?
  inline bool is_hyperon() const {
    const auto abs_code = std::abs(code());
    switch (abs_code) {
      case pdg::Lambda:
      case pdg::Sigma_p:
      case pdg::Sigma_z:
      case pdg::Sigma_m:
      case pdg::Xi_z:
      case pdg::Xi_m:
      case pdg::Omega_m:
        return true;
      default:
        return false;
    }
  }
  /// Is this a Xi(1321)?
  inline bool is_xi1321() const {
    const auto abs_code = std::abs(code());
    return (abs_code == pdg::Xi_z) || (abs_code == pdg::Xi_m);
  }
  /// Is this a Omega(1672)?
  inline bool is_Omega1672() const {
    const auto abs_code = std::abs(code());
    return (abs_code == pdg::Omega_m);
  }
  /// Is this a kaon (K+, K-, K0, Kbar0)?
  inline bool is_kaon() const {
    const auto abs_code = std::abs(code());
    return (abs_code == pdg::K_p) || (abs_code == pdg::K_z);
  }
  /// Is this a pion (pi+/pi0/pi-)?
  inline bool is_pion() const {
    const auto c = code();
    return (c == pdg::pi_z) || (c == pdg::pi_p) || (c == pdg::pi_m);
  }
  /// Is this a rho meson (rho+/rho0/rho-)?
  inline bool is_rho() const {
    const auto c = code();
    return (c == pdg::rho_z) || (c == pdg::rho_p) || (c == pdg::rho_m);
  }

  /** Determine whether a particle has a distinct antiparticle
    * (or whether it is its own antiparticle). */
  bool has_antiparticle() const {
    if (is_hadron()) {
      return (baryon_number() != 0) || (digits_.n_q2_ != digits_.n_q3_);
    } else {
      return digits_.n_q3_ == 1;  // leptons!
    }
  }
  /** returns twice the isospin-3 component \f$I_3\f$.
   *
   * This is calculated from the sum of net_quark_number of up and down.
   */
  inline int isospin3() const {
    // net_quark_number(2) is the number of u quarks,
    // net_quark_number(1) is the number of d quarks.
    return net_quark_number(2) - net_quark_number(1);
  }
  /** returns the net number of \f$\bar s\f$ quarks.
   *
   * For particles with one strange quark, -1 is returned.
   */
  inline int strangeness() const { return -net_quark_number(3); }
  /** returns the net number of \f$c\f$ quarks
   *
   * For particles with one charm quark, +1 is returned.
   */
  inline int charmness() const { return +net_quark_number(4); }
  /** returns the net number of \f$\bar b\f$ quarks
   *
   * For particles with one bottom quark, -1 is returned.
   */
  inline int bottomness() const { return -net_quark_number(5); }
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
      /* This loops over d,u,s,c,b,t quarks (the latter can be safely ignored,
       * but I don't think this will be a bottle neck. */
      for (int i = 1; i < 7; i++) {
        /* u,c,t quarks have charge = 2/3 e, while d,s,b quarks have -1/3 e.
         * The antiparticle sign is already in net_quark_number. */
        Q += (i % 2 == 0 ? 2 : -1) * net_quark_number(i);
      }
      return Q / 3;
    }
    /* non-hadron:
     * Leptons: 11, 13, 15 are e, μ, τ and have a charge -1, while
     *          12, 14, 16 are the neutrinos that have no charge. */
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
      if (digits_.n_J_ == 0) {
        return 0;  // special cases: K0_L=0x130 & K0_S=0x310
      } else {
        return digits_.n_J_ - 1;
      }
    }
    // this assumes that we only have white particles (no single
    // quarks): Electroweak fermions have 11-17, so the
    // second-to-last-digit is the spin. The same for the Bosons: they
    // have 21-29 and 2spin = 2 (this fails for the Higgs).
    return digits_.n_q3_;
  }
  /** Returns the spin degeneracy \f$2s + 1\f$ of a particle **/
  inline unsigned int spin_degeneracy() const {
    if (is_hadron() && digits_.n_J_ > 0) {
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
    if (!is_hadron()) {
      return 0;
    }
    return chunks_.quarks_;
  }

  /** Returns quark content as an array, useful for interfacing to Pythia.
   *  The return is always an array of three numbers, which are pdgcodes
   *  of quarks: 1 - d, 2 - u, 3 - s, 4 - c, 5 - b. Antiquarks get a negative
   *  sign. For mesons the first number in array is always 0.
   *
   *  There is a difficulty with mesons that are a superposition, for example
   *  \f$ \pi^0 = \frac{1}{\sqrt{2}}(u \bar{u} + d \bar{d}) \f$. Currently for
   *  \f$ \pi^0 \f$ just {0, 1, -1} is returned.
   *
   */
  std::array<int, 3> quark_content() const {
    std::array<int, 3> result = { static_cast<int>(digits_.n_q1_),
                                  static_cast<int>(digits_.n_q2_),
                                  static_cast<int>(digits_.n_q3_)};
    if (is_hadron()) {
      // Antibaryons
      if (digits_.n_q1_ != 0 && digits_.antiparticle_) {
        for (size_t i = 0; i < 3; i++) {
          result[i] = -result[i];
        }
      }
      // Mesons
      if (digits_.n_q1_ == 0) {
        // Own antiparticle
        if (digits_.n_q2_ == digits_.n_q3_) {
          result[2] = -result[2];
        } else {
          // Like pi-
          if (digits_.antiparticle_) {
            result[1] = -result[1];
          // Like pi+
          } else {
            result[2] = -result[2];
          }
        }
      }
    } else {
      result = {0, 0, 0};
    }
    return result;
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
    return dump_ < rhs.dump_;
    // the complex thing to do here is to calculate:
    //   code() < rhs.code()
    // but for getting a total order that's overkill. The uint32_t value in
    // dump_ works just fine.
  }
  /// returns if the codes are equal
  inline bool operator==(const PdgCode rhs) const { return dump_ == rhs.dump_; }
  /// returns if the codes are not equal.
  inline bool operator!=(const PdgCode rhs) const { return !(*this == rhs); }
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

  /** returns an integer with decimal representation of the code.
   *
   * This is necessary for ROOT output.
   *
   */
  int get_decimal() const {
    return antiparticle_sign() *
           (digits_.n_J_ + digits_.n_q3_ * 10 + digits_.n_q2_ * 100 +
            digits_.n_q1_ * 1000 + digits_.n_L_ * 10000 +
            digits_.n_R_ * 100000 + digits_.n_ * 1000000);
  }

 private:
// amend this line with something that identifies your compiler if its
// bit field order is like in the gnu c compiler for 64 bit
// architectures (if you are unsure, try one and check the pdgcode
// test).
#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__)) || \
    defined(DOXYGEN)
#define SMASH_BITFIELD_ORDER_ 1
// put your compiler here if the bit field order is reversed w.r.t. gnu
// c compiler for 64 bit.
#elif defined(__OTHER_COMPILER__)
#define SMASH_BITFIELD_ORDER_ 2
#else
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
#if SMASH_BITFIELD_ORDER_ == 1
      /// spin quantum number \f$n_J = 2 J + 1\f$.
      std::uint32_t n_J_ : 4;
      /// third quark field
      std::uint32_t n_q3_ : 4;
      /// second quark field
      std::uint32_t n_q2_ : 4;
      /// first quark field. 0 for mesons.
      std::uint32_t n_q1_ : 4;
      /// "angular momentum"
      std::uint32_t n_L_ : 4;
      /// "radial excitation"
      std::uint32_t n_R_ : 4;
      /// first field: "counter"
      std::uint32_t n_ : 4, : 3;
      /// first bit: stores the sign.
      bool antiparticle_ : 1;
#else  // reverse ordering
      bool antiparticle_ : 1, : 3;
      std::uint32_t n_ : 4;
      std::uint32_t n_R_ : 4;
      std::uint32_t n_L_ : 4;
      std::uint32_t n_q1_ : 4;
      std::uint32_t n_q2_ : 4;
      std::uint32_t n_q3_ : 4;
      std::uint32_t n_J_ : 4;
#endif
    } digits_;
    /** the bitfield dumped into a single integer. Please note that the
     * 2nd, 3rd and 4th highest bits are possibly undefined.
     **/
    std::uint32_t dump_;
    /** chunk collection: here, the chunks with \f$nn_Rn_L\f$ and
     * \f$n_{q_1}n_{q_2}n_{q_3}\f$ are directly accessible.
     */
    struct {
#if SMASH_BITFIELD_ORDER_ == 1
      std::uint32_t : 4;
      /// the quark digits n_q{1,2,3}_
      std::uint32_t quarks_ : 12;
      /// the excitation digits n_, n_R_, n_L_
      std::uint32_t excitation_ : 12, : 4;
#else  // reverse ordering
      std::uint32_t : 4, excitation_ : 12;
      std::uint32_t quarks_ : 12, : 4;
#endif
    } chunks_;
  };

  /** Returns an unsigned integer with the PDG code in hexadecimal
   *  (disregarding the antiparticle flag). */
  inline std::uint32_t ucode() const { return (dump_ & 0x0fffffff); }
  /** returns the net number of quarks with given flavour number
   *
   * \param quark PDG Code of quark: (1..6) = (d,u,s,c,b,t)
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
      throw InvalidPdgCode("PdgCode: Invalid character " +
                           std::string(&inp, 1) + " found.\n");
    }
    // the last four digits are the number; they should not be > 9
    // (i.e., one of [:;<=>?])
    if ((inp & 0x0f) > 9) {
      throw InvalidPdgCode("PdgCode: Invalid digit " + std::string(&inp, 1) +
                           " found.\n");
    }
    // now that we've checked that the first bits are correct and the
    // last bits are a number, we can return the last bits.
    return (inp & 0x0f);
  }

  /// takes a string and sets the fields.
  inline void set_from_string(const std::string& codestring) {
    dump_ = 0;
    // implicit with the above: digits_.antiparticle_ = false;
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
      if (digits_.n_q1_ > 6) {
        throw InvalidPdgCode("Invalid PDG code " + codestring + " (n_q1>6)");
      }
    }
    // 3rd from last is n_q2_.
    if (length > 2 + sign) {
      digits_.n_q2_ = get_digit_from_char(codestring[c++]);
      if (digits_.n_q2_ > 6) {
        throw InvalidPdgCode("Invalid PDG code " + codestring + " (n_q2>6)");
      }
    }
    // next to last is n_q3_.
    if (length > 1 + sign) {
      digits_.n_q3_ = get_digit_from_char(codestring[c++]);
      if (digits_.n_q3_ > 6) {
        throw InvalidPdgCode("Invalid PDG code " + codestring + " (n_q3>6)");
      }
    }
    // last digit is the spin degeneracy.
    if (length > sign) {
      digits_.n_J_ = get_digit_from_char(codestring[c++]);
    } else {
      throw InvalidPdgCode(
          "String \"" + codestring +
          "\" only consists of a sign, that is no valid PDG Code\n");
    }
    check();
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
                           " in PDG Code " + string());
    }
    check();
  }
};

static_assert(sizeof(PdgCode) == 4, "should fit into 32 bit integer");

std::istream& operator>>(std::istream& is, PdgCode& code);
/**\ingroup logging
 * Writes the textual representation of the PDG code to the output stream.
 */
std::ostream& operator<<(std::ostream& is, const PdgCode& code);

/** Checks if two given particles represent a lepton pair (e+e- or mu+mu-). */
inline bool is_dilepton(const PdgCode pdg1, const PdgCode pdg2) {
  const auto c1 = pdg1.code();
  const auto c2 = pdg2.code();
  const auto min = std::min(c1, c2);
  const auto max = std::max(c1, c2);
  return (max == 0x11 && min == -0x11) || (max == 0x13 && min == -0x13);
}

/** Checks if two of the three given particles represent a lepton pair
 * (e+e- or mu+mu-).*/
inline bool has_lepton_pair(const PdgCode pdg1, const PdgCode pdg2,
                            const PdgCode pdg3) {
  return is_dilepton(pdg1, pdg2) || is_dilepton(pdg1, pdg3) ||
         is_dilepton(pdg2, pdg3);
}

}  // namespace Smash

#endif  // SRC_INCLUDE_PDGCODE_H_
