/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_PDGCODE_H_
#define SRC_INCLUDE_PDGCODE_H_

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <iosfwd>
#include <sstream>
#include <stdexcept>
#include <string>

#include <iostream>

#include "pdgcode_constants.h"

namespace smash {

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
 * Representing nuclei
 * -------------------
 *
 * Following PDG standard, nuclei are represented by codes ±10LZZZAAAI, where
 * L is number of Lambdas inside the nucleus, ZZZ is charge, AAA is mass
 * number and I is used for excitations. Internally nuclei are represented
 * in a different way from hadrons, but all accessors (charge, baryon number,
 * etc) work in the same way.
 *
 * Normally nuclei in SMASH are simulated as a collection of protons and
 * neutrons, so there is no need in their PDG codes. However, it is
 * interesting to study light nuclei production, considering them as single
 * pointlike hadrons. This justifies introduction of nuclear PDG codes here.
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
 * Also, tetra- and pentaquarks cannot be represented; that, though,
 * is a problem of the PDG Numbering Scheme rather than of this class.
 */

class PdgCode {
 public:
  /**
   * \ingroup exception
   * thrown for invalid inputs
   */
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
  /**
   * Initialize using a string
   * The string is interpreted as a hexadecimal number, i.e., \c 211 is
   * interpreted as \c 0x211 = \f$529_{10}\f$.
   */
  explicit PdgCode(const std::string& codestring) {
    set_from_string(codestring);
  }

  /**
   * Receive a signed integer and process it into a PDG Code. The sign
   * is taken as antiparticle boolean, while the absolute value of the
   * integer is used as hexdigits.
   * \param[in] codenumber a signed integer which represent the PDG code
   * The number 0x221 is interpreted as an η meson,
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
  /**
   * receive an unsigned integer and process it into a PDG Code. The
   * first bit is taken and used as antiparticle boolean.
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

  /**
   * Checks the integer for invalid hex digits.
   *
   * Usually all digits are at least <= 9. The n_q digits are even <= 6
   * (because there are only six quarks). The only exception is n_J, where
   * we allow f = 15, which is the largest hexadecimal digit.
   * If one of the hex digits is not also a valid decimal digit,
   * something possibly went wrong - maybe some user of this class forgot
   * to prefix the input with '0x' and thus passed 221 instead of 0x221.
   * \return a bitmask indicating the offending digits. In the above
   * example, 221 = 0xd3, the second-to-last-digit is the offending one,
   * to the return value is 0b10 = 0x2 = 2.
   */
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
    if (digits_.n_J_ > 15) {
      fail |= 1;
    }
    return fail;
  }

  /**
   * Do all sorts of validity checks.
   * \throw InvalidPdgCode if meson has even n_J_ (fermionic spin)
   * \throw InvalidPdgCode if baryon has odd n_J_ (bosonic spin)
   * \throw InvalidPdgCode if n_J_ is 0 (spin is not defined.)
   * \throw InvalidPdgCode if particle does not have antiparticle when
   * it is supposed to do.
   */
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

  /// Dumps the bitfield into an unsigned integer.
  inline std::uint32_t dump() const {
    // this cuts the three unused bits.
    return (dump_ & 0x8fffffff);
  }

  /// \return a signed integer with the PDG code in hexadecimal.
  inline std::int32_t code() const { return antiparticle_sign() * ucode(); }

  /// \return the PDG Code as a decimal string.
  inline std::string string() const {
    std::stringstream ss;
    ss << get_decimal();
    return ss.str();
  }

  /// Construct the antiparticle to a given PDG code.
  PdgCode get_antiparticle() const {
    PdgCode result = *this;
    result.digits_.antiparticle_ = !digits_.antiparticle_;
    return result;
  }

  /**
   * Construct PDG code from decimal number.
   * \param[in] pdgcode_decimal decimal integer representing the PDG code
   */
  static PdgCode from_decimal(const int pdgcode_decimal) {
    // Nucleus
    if (std::abs(pdgcode_decimal) > 1E9) {
      return PdgCode(std::to_string(pdgcode_decimal));
    }
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

  /// \return true if this is a nucleus, false otherwise
  inline bool is_nucleus() const {
    assert(digits_.is_nucleus_ == nucleus_.is_nucleus_);
    return nucleus_.is_nucleus_;
  }

  /// \return true if this is a baryon, antibaryon or meson.
  inline bool is_hadron() const {
    return (digits_.n_q3_ != 0 && digits_.n_q2_ != 0 && !is_nucleus());
  }

  /// \return true if this is a lepton.
  inline bool is_lepton() const {
    return (digits_.n_q1_ == 0 && digits_.n_q2_ == 0 && digits_.n_q3_ == 1 &&
            !is_nucleus());
  }

  /// \return the baryon number of the particle.
  inline int baryon_number() const {
    if (is_nucleus()) {
      return static_cast<int>(nucleus_.A_) * antiparticle_sign();
    }
    if (!is_hadron() || digits_.n_q1_ == 0) {
      return 0;
    }
    return antiparticle_sign();
  }
  /// \return whether this PDG code identifies a baryon.
  inline bool is_baryon() const { return is_hadron() && digits_.n_q1_ != 0; }

  /// \return whether this PDG code identifies a meson.
  inline bool is_meson() const { return is_hadron() && digits_.n_q1_ == 0; }

  /// \return whether this is a nucleon/anti-nucleon (p, n, -p, -n)
  inline bool is_nucleon() const {
    const auto abs_code = std::abs(code());
    return (abs_code == pdg::p || abs_code == pdg::n);
  }

  /// \return whether this is a N*(1535) (+/0)
  inline bool is_Nstar1535() const {
    const auto abs_code = std::abs(code());
    return (abs_code == pdg::N1535_p || abs_code == pdg::N1535_z);
  }

  /// \return whether this is a Delta(1232) (with anti-delta)
  inline bool is_Delta() const {
    const auto abs_code = std::abs(code());
    return (abs_code == pdg::Delta_pp || abs_code == pdg::Delta_p ||
            abs_code == pdg::Delta_z || abs_code == pdg::Delta_m);
  }

  /// \return whether this is a hyperon (Lambda, Sigma, Xi, Omega)
  inline bool is_hyperon() const { return is_hadron() && digits_.n_q1_ == 3; }

  /// \return whether this is a Omega baryon
  inline bool is_Omega() const {
    return is_hyperon() && digits_.n_q2_ == 3 && digits_.n_q3_ == 3;
  }

  /// \return whether this is a Xi baryon
  inline bool is_Xi() const {
    return is_hyperon() && digits_.n_q2_ == 3 && digits_.n_q3_ != 3;
  }

  /// \return whether this is a Lambda baryon
  inline bool is_Lambda() const {
    return is_hyperon() && digits_.n_q2_ == 1 && digits_.n_q3_ == 2;
  }

  /// \return whether this is a Sigma baryon
  inline bool is_Sigma() const {
    return is_hyperon() && digits_.n_q2_ != 3 && !is_Lambda();
  }

  /// \return whether this is a kaon (K+, K-, K0, Kbar0)
  inline bool is_kaon() const {
    const auto abs_code = std::abs(code());
    return (abs_code == pdg::K_p) || (abs_code == pdg::K_z);
  }

  /// \return whether this is a pion (pi+/pi0/pi-)
  inline bool is_pion() const {
    const auto c = code();
    return (c == pdg::pi_z) || (c == pdg::pi_p) || (c == pdg::pi_m);
  }

  /// \return whether this is a rho meson (rho+/rho0/rho-)
  inline bool is_rho() const {
    const auto c = code();
    return (c == pdg::rho_z) || (c == pdg::rho_p) || (c == pdg::rho_m);
  }

  /// \return whether this is (anti-)deuteron
  inline bool is_deuteron() const {
    const int dec = get_decimal();
    return is_nucleus() && (dec == pdg::decimal_d || dec == pdg::decimal_antid);
  }

  /**
   * \return whether a particle has a distinct antiparticle
   * (or whether it is its own antiparticle).
   */
  bool has_antiparticle() const {
    if (is_nucleus()) {
      return true;
    }
    if (is_hadron()) {
      return (baryon_number() != 0) || (digits_.n_q2_ != digits_.n_q3_);
    } else {
      return digits_.n_q3_ == 1;  // leptons!
    }
  }

  /**
   * \return twice the isospin-3 component \f$I_3\f$.
   *
   * This is calculated from the sum of net_quark_number of up and down.
   */
  inline int isospin3() const {
    /* net_quark_number(2) is the number of u quarks,
     * net_quark_number(1) is the number of d quarks. */
    return net_quark_number(2) - net_quark_number(1);
  }

  /**
   * \return the fraction number of strange quarks
   *         (strange + anti-strange) / total
   *
   * This is useful for the AQM cross-section scaling, and needs to
   * be positive definite.
   */
  inline double frac_strange() const {
    /* The quarkonium state has 0 net strangeness
     *  but there are actually 2 strange quarks out of 2 total */
    if (is_hadron() && digits_.n_q3_ == 3 && digits_.n_q2_ == 3) {
      return 1.;
    } else {
      // For all other cases, there isn't both a strange and anti-strange
      if (is_baryon()) {
        return abs(strangeness()) / 3.;
      } else if (is_meson()) {
        return abs(strangeness()) / 2.;
      } else {
        /* If not baryon or meson, this should be 0, as AQM does not
         * extend to non-hadrons */
        return 0.;
      }
    }
  }

  /**
   * \return the net number of \f$\bar s\f$ quarks.
   *
   * For particles with one strange quark, -1 is returned.
   */
  inline int strangeness() const { return -net_quark_number(3); }

  /**
   * \return the net number of \f$c\f$ quarks
   *
   * For particles with one charm quark, +1 is returned.
   */
  inline int charmness() const { return +net_quark_number(4); }

  /**
   * \return the net number of \f$\bar b\f$ quarks
   *
   * For particles with one bottom quark, -1 is returned.
   */
  inline int bottomness() const { return -net_quark_number(5); }

  /**
   * The charge of the particle.
   * The charge is calculated from the quark content (for hadrons) or
   * basically tabulated; currently leptons, neutrinos and the standard
   * model gauge bosons are known; unknown particles return a charge of
   * 0.
   * \return charge of the particle
   */
  int charge() const {
    if (is_hadron() || is_nucleus()) {
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
    /* Bosons: 24 is the W+, all else is uncharged.
     * we ignore the first digits so that this also finds strange gauge
     * boson "resonances" (in particular, \f$\tilde \chi_1^+\f$ with PDG
     * Code 1000024). */
    if ((dump_ & 0x0000ffff) == 0x24) {
      return antiparticle_sign();
    }
    // default (this includes all other Bosons) is 0.
    return 0;
  }

  /**
   * \todo (oliiny): take care of spin for nuclei
   * \return twice the spin of a particle.
   *
   * The code is good for hadrons, leptons and spin-1-bosons. It returns
   * 2 (meaning spin=1) for the Higgs, though.
   */
  inline unsigned int spin() const {
    if (is_nucleus()) {
      /* Currently the only nucleus I care about is deutron,
       * which has spin one. */
      return 2;
    }

    if (is_hadron()) {
      if (digits_.n_J_ == 0) {
        return 0;  // special cases: K0_L=0x130 & K0_S=0x310
      } else {
        return digits_.n_J_ - 1;
      }
    }
    /* this assumes that we only have white particles (no single
     * quarks): Electroweak fermions have 11-17, so the
     * second-to-last-digit is the spin. The same for the Bosons: they
     * have 21-29 and 2spin = 2 (this fails for the Higgs). */
    return digits_.n_q3_;
  }
  /// \return the spin degeneracy \f$2s + 1\f$ of a particle.
  inline unsigned int spin_degeneracy() const {
    if (is_hadron() && digits_.n_J_ > 0) {
      return digits_.n_J_;
    }
    return spin() + 1;
  }
  /// \return -1 for antiparticles and +1 for particles.
  inline int antiparticle_sign() const {
    return (digits_.antiparticle_ ? -1 : +1);
  }
  /// \return an integer with only the quark numbers set.
  inline std::int32_t quarks() const {
    if (!is_hadron() || is_nucleus()) {
      return 0;
    }
    return chunks_.quarks_;
  }

  /**
   * The return is always an array of three numbers, which are pdgcodes
   * of quarks: 1 - d, 2 - u, 3 - s, 4 - c, 5 - b. Antiquarks get a negative
   * sign. For mesons the first number in array is always 0.
   * There is a difficulty with mesons that are a superposition, for example
   * \f$ \pi^0 = \frac{1}{\sqrt{2}}(u \bar{u} + d \bar{d}) \f$. Currently for
   * \f$ \pi^0 \f$ just {0, 1, -1} is returned.
   * \return quark content as an array.
   */
  std::array<int, 3> quark_content() const {
    std::array<int, 3> result = {static_cast<int>(digits_.n_q1_),
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
        // add extra minus sign according to the pdg convention
        if (digits_.n_q2_ != digits_.n_q3_ && digits_.n_q2_ % 2 == 1) {
          for (int i = 1; i <= 2; i++) {
            result[i] = -result[i];
          }
        }
      }
    } else {
      result = {0, 0, 0};
    }
    return result;
  }

  /**
   * \return whether a particle contains at least the given number of
   * valence quarks.
   * \param[in] valence_quarks_required number of valence quarks
   * that particle is supposed to contain.
   *
   * \throw std::runtime_error if it is not a hadron
   *
   * This is necessary for string fragmentation.
   */
  bool contains_enough_valence_quarks(int valence_quarks_required) const;

  /****************************************************************************
   *                                                                          *
   * operations with more than one PDG Code                                   *
   *                                                                          *
   ****************************************************************************/

  /**
   * Sorts PDG Codes according to their numeric value.
   * This is used by std::map
   */
  inline bool operator<(const PdgCode rhs) const {
    return dump_ < rhs.dump_;
    /* the complex thing to do here is to calculate:
     *   code() < rhs.code()
     * but for getting a total order that's overkill. The uint32_t value in
     * dump_ works just fine. */
  }

  /// \return if the codes are equal
  inline bool operator==(const PdgCode rhs) const { return dump_ == rhs.dump_; }

  /// \return if the codes are not equal.
  inline bool operator!=(const PdgCode rhs) const { return !(*this == rhs); }

  /// \return if the code of rhs is the inverse of this one.
  inline bool is_antiparticle_of(const PdgCode rhs) const {
    return code() == -rhs.code();
  }

  /// istream >> PdgCode assigns the PDG Code from an istream.
  friend std::istream& operator>>(std::istream& is, PdgCode& code);

  /**
   * PdgCode 0x0 is guaranteed not to be valid by the PDG standard, but
   * it passes all tests here, so we can use it to show some code is not
   * yet set.
   */
  static PdgCode invalid() { return PdgCode(0x0); }

  /**
   * \return an integer with decimal representation of the code.
   * If the spin is too large for the last digit, an additional digit at the
   * beginning will be used, so that the sum of the first and the last digit is
   * the spin.
   * This is used for binary and ROOT output.
   *
   * \throw InvalidPdgCode if the spin degeneracy is larger than 9
   */
  int get_decimal() const {
    if (is_nucleus()) {
      // ±10LZZZAAAI
      return antiparticle_sign() *
             (nucleus_.I_ + 10 * nucleus_.A_ + 10000 * nucleus_.Z_ +
              10000000 * nucleus_.n_Lambda_ + 1000000000);
    }
    int n_J_1 = 0;
    int n_J_2 = digits_.n_J_;
    if (n_J_2 > 9) {
      n_J_1 = n_J_2 - 9;
      n_J_2 = 9;
      if (n_J_2 > 9) {
        throw InvalidPdgCode("n_J is too large\n");
      }
    }
    return antiparticle_sign() *
           (n_J_2 + digits_.n_q3_ * 10 + digits_.n_q2_ * 100 +
            digits_.n_q1_ * 1000 + digits_.n_L_ * 10000 +
            digits_.n_R_ * 100000 + digits_.n_ * 1000000 + n_J_1 * 10000000);
  }

  /// Remove all excitation, except spin. Sign and quark content remains.
  void deexcite() {
    if (!is_nucleus()) {
      chunks_.excitation_ = 0;
    } else {
      nucleus_.I_ = 0;
    }
  }

  /**
   * Returns the net number of quarks with given flavour number
   * For public use, see strangeness(), charmness(), bottomness() and
   * isospin3().
   * \param[in] quark PDG Code of quark: (1..6) = (d,u,s,c,b,t)
   * \return for the net number of quarks (\#quarks - \#antiquarks)
   *
   * \throw std::invalid_argument
   * if quark is not any of d, u, s, c, b and t quarks
   */
  int net_quark_number(const int quark) const;

 private:
/* amend this line with something that identifies your compiler if its
 * bit field order is like in the gnu c compiler for 64 bit
 * architectures (if you are unsure, try one and check the pdgcode
 * test). */
#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__)) || \
    defined(DOXYGEN)
#define SMASH_BITFIELD_ORDER_ 1
/* put your compiler here if the bit field order is reversed w.r.t. gnu
 * c compiler for 64 bit. */
#elif defined(__OTHER_COMPILER__)
#define SMASH_BITFIELD_ORDER_ 2
#else
#error "Please determine the correct bit-field order for your target/compiler"
#endif
  /**
   * The union holds the data; either as a single integer dump_, as a
   * single-digit bitfield digits_ or as a multiple-digits bitfield
   * chunks_.
   */
  union {
    /**
     * The single digits collection of the code. Here, every PDG code
     * digits is directly accessible.
     */
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
      std::uint32_t n_ : 4, : 2;
      /// 1 for nuclei, 0 for the rest
      bool is_nucleus_ : 1;
      /// first bit: stores the sign.
      bool antiparticle_ : 1;
#else  // reverse ordering
      bool antiparticle_ : 1;
      bool is_nucleus : 1, : 2;
      std::uint32_t n_ : 4;
      std::uint32_t n_R_ : 4;
      std::uint32_t n_L_ : 4;
      std::uint32_t n_q1_ : 4;
      std::uint32_t n_q2_ : 4;
      std::uint32_t n_q3_ : 4;
      std::uint32_t n_J_ : 4;
#endif
    } digits_;
    /**
     * The bitfield dumped into a single integer. Please note that the
     * 2nd, 3rd and 4th highest bits are possibly undefined.
     */
    std::uint32_t dump_;
    /**
     * Chunk collection: here, the chunks with \f$nn_Rn_L\f$ and
     * \f$n_{q_1}n_{q_2}n_{q_3}\f$ are directly accessible.
     */
    struct {
#if SMASH_BITFIELD_ORDER_ == 1
      std::uint32_t : 4;
      /// The quark digits n_q{1,2,3}_
      std::uint32_t quarks_ : 12;
      /// The excitation digits n_, n_R_, n_L_
      std::uint32_t excitation_ : 12, : 4;
#else  /// Reverse ordering
      std::uint32_t : 4, excitation_ : 12;
      std::uint32_t quarks_ : 12, : 4;
#endif
    } chunks_;
    /// Structure for the nuclei
    struct {
#if SMASH_BITFIELD_ORDER_ == 1
      std::uint32_t n_Lambda_ : 6;
      std::uint32_t Z_ : 10;
      std::uint32_t A_ : 10;
      std::uint32_t I_ : 4;
      bool is_nucleus_ : 1;
      bool antiparticle_ : 1;
#else  // Reverse ordering
      bool antiparticle_ : 1;
      bool is_nucleus_ : 1;
      std::uint32_t I_ : 4;
      std::uint32_t A_ : 10;
      std::uint32_t Z_ : 10;
      std::uint32_t n_Lambda_ : 6;
#endif
    } nucleus_;
  };

  /**
   * \return an unsigned integer with the PDG code in hexadecimal
   * (disregarding the antiparticle flag).
   */
  inline std::uint32_t ucode() const { return (dump_ & 0x0fffffff); }

  /**
   * \return digits from a hexadecimal character.
   * \param[in] inp character which is translated into digit
   *
   * \throw InvalidPdgCode if character does not correspond to digit
   */
  inline std::uint32_t get_digit_from_char(const char inp) const {
    // Decimal digit
    if (48 <= inp && inp <= 57) {
      return inp - 48;
    }
    // Hexdecimal digit, uppercase
    if (65 <= inp && inp <= 70) {
      return inp - 65 + 10;
    }
    // Hexdecimal digit, lowercase
    if (97 <= inp && inp <= 102) {
      return inp - 97 + 10;
    }
    throw InvalidPdgCode("PdgCode: Invalid character " + std::string(&inp, 1) +
                         " found.\n");
  }

  /**
   * Set the PDG code from the given string.
   * This supports hexdecimal digits. If the last digit is not enough to
   * represent the spin, a digit can be added at the beginning which will be
   * added to the total spin.
   * \param[in] codestring string which is translated into PdgCode
   *
   * \throw InvalidPdgCode if the input string is empty
   * \throw InvalidPdgCode
   * if it is a nucleus whose PDG code does not begin with 10
   * \throw InvalidPdgCode
   * if it is not a nucleus while number of digits is more than 8
   * \throw InvalidPdgCode
   * if the 1st quark field is not any of d, u, s, c, b and t quarks
   * \throw InvalidPdgCode
   * if the 2nd quark field is not any of d, u, s, c, b and t quarks
   * \throw InvalidPdgCode
   * if the 3rd quark field is not any of d, u, s, c, b and t quarks
   * \throw InvalidPdgCode
   * if there is nothing else but sign
   */
  inline void set_from_string(const std::string& codestring) {
    dump_ = 0;
    // Implicit with the above: digits_.antiparticle_ = false;
    digits_.n_ = digits_.n_R_ = digits_.n_L_ = digits_.n_q1_ = digits_.n_q2_ =
        digits_.n_q3_ = digits_.n_J_ = digits_.is_nucleus_ = 0;
    size_t length = codestring.size();
    if (length < 1) {
      throw InvalidPdgCode("Empty string does not contain PDG Code\n");
    }
    int c = 0;
    /* Look at current character; if it is a + or minus sign, read it
     * and advance to next char. */
    if (codestring[c] == '-') {
      digits_.antiparticle_ = true;
      ++c;
    } else if (codestring[c] == '+') {
      digits_.antiparticle_ = false;
      ++c;
    }
    // Save if the first character was a sign:
    unsigned int sign = c;

    // Nucleus
    if (length == 10 + sign) {
      nucleus_.is_nucleus_ = true;
      if (codestring[c] != '1' || codestring[c + 1] != '0') {
        throw InvalidPdgCode("Pdg code of nucleus \"" + codestring +
                             "\" should start with 10\n");
      }
      c += 2;
      // ±10LZZZAAAI is the standard for nuclei
      std::array<int, 8> digits;
      for (int i = 0; i < 8; i++) {
        digits[i] = get_digit_from_char(codestring[c + i]);
      }
      nucleus_.n_Lambda_ = digits[0];
      nucleus_.Z_ = 100 * digits[1] + 10 * digits[2] + digits[3];
      nucleus_.A_ = 100 * digits[4] + 10 * digits[5] + digits[6];
      nucleus_.I_ = digits[7];
      return;
    }

    // Codestring shouldn't be longer than 8 + sign, except for nuclei
    if (length > 8 + sign) {
      throw InvalidPdgCode("String \"" + codestring +
                           "\" too long for PDG Code\n");
    }
    /* Please note that in what follows, we actually need c++, not ++c.
     * first digit is used for n_J if the last digit is not enough. */
    if (length > 7 + sign) {
      digits_.n_J_ += get_digit_from_char(codestring[c++]);
    }
    // Codestring has 7 digits? 7th from last goes in n_.
    if (length > 6 + sign) {
      digits_.n_ = get_digit_from_char(codestring[c++]);
    }
    // It has 6 or 7 digits? 6th from last is n_R_.
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
    // Next to last is n_q3_.
    if (length > 1 + sign) {
      digits_.n_q3_ = get_digit_from_char(codestring[c++]);
      if (digits_.n_q3_ > 6) {
        throw InvalidPdgCode("Invalid PDG code " + codestring + " (n_q3>6)");
      }
    }
    // Last digit is the spin degeneracy.
    if (length > sign) {
      digits_.n_J_ += get_digit_from_char(codestring[c++]);
    } else {
      throw InvalidPdgCode(
          "String \"" + codestring +
          "\" only consists of a sign, that is no valid PDG Code\n");
    }
    check();
  }

  /**
   * Sets the bitfield from an unsigned integer. Usually called from
   * the constructors.
   * \param[in] abscode integer which replace PDG code except sign
   *
   * \throw InvalidPdgCode if input is not a valid PDG code
   *
   * \see PdgCode::test_code
   */
  inline void set_fields(std::uint32_t abscode) {
    /* "dump_ =" overwrites antiparticle_, but this needs to have been set
     * already, so we carry it around the assignment. */
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

/**
 * Sets the PDG code from the textual representation
 * in the input stream.
 * \param[in] is input string
 * \param[out] code PdgCode to be set
 */
std::istream& operator>>(std::istream& is, PdgCode& code);
/**
 * \ingroup logging
 * Writes the textual representation of the PDG code
 * to the output stream.
 */
std::ostream& operator<<(std::ostream& is, const PdgCode& code);

/// \return if two given particles represent a lepton pair (e+e- or mu+mu-).
inline bool is_dilepton(const PdgCode pdg1, const PdgCode pdg2) {
  const auto c1 = pdg1.code();
  const auto c2 = pdg2.code();
  const auto min = std::min(c1, c2);
  const auto max = std::max(c1, c2);
  return (max == 0x11 && min == -0x11) || (max == 0x13 && min == -0x13);
}

/**
 * \return if two of the three given particles represent a lepton pair
 * (e+e- or mu+mu-).
 */
inline bool has_lepton_pair(const PdgCode pdg1, const PdgCode pdg2,
                            const PdgCode pdg3) {
  return is_dilepton(pdg1, pdg2) || is_dilepton(pdg1, pdg3) ||
         is_dilepton(pdg2, pdg3);
}

}  // namespace smash

#endif  // SRC_INCLUDE_PDGCODE_H_
