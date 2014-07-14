/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_QUANTUMNUMBERS_H_
#define SRC_INCLUDE_QUANTUMNUMBERS_H_

#include<iostream>
#include<sstream>

#include "fourvector.h"
#include "numerics.h"
#include "particles.h"
#include "pdgcode.h"

namespace Smash {

/** QuantumNumbers - a container for storing conserved values
 *
 * QuantumNumbers can be used to compare all values that should be
 * conserved during the evolution, in particular for the whole system at
 * fixed intervals, or for a subset of particles before- and after a
 * DecayAction / ScatterAction / Action.
 *
 * QuantumNumbers can also be used to store, and retrieve, the total
 * value of the quantum numbers in e.g. a ScatterAction, before
 * distributing it to the new particles.
 *
 * Currently, momentum conservation (including energy conservation),
 * charge-, isospin3-, (net-) strangeness-, (net-) bottomness-, (net-)
 * charmness- and (net-) baryon number conservation are checked. (It
 * should be noted, or repeated, that only the net quantities are
 * conserved, hence this is what is stored and compared with this
 * class.)
 *
 * Usage
 * -----
 *
 * \code
 * QuantumNumbers before(particlelist);
 * do_something_with(particlelist);
 * QuantumNumbers after(particlelist);
 * printf("%s", before.report_deviations(after).c_str());
 * if (before != after) {
 *   throw std::runtime_error(before.report_deviations(after));
 * }
 * \endcode
 */
class QuantumNumbers {
 public:
  /// construct QuantumNumbers collection with all fields 0.
  QuantumNumbers()
        : momentum_(0,0,0,0),
          charge_(0),
          isospin3_(0),
          strangeness_(0),
          charmness_(0),
          bottomness_(0),
          baryon_number_(0) {}
  /** construct QuantumNumbers collection
   *
   * \param m Momentum FourVector
   * \param q Charge
   * \param i3 Isospin
   * \param s Strangeness
   * \param c Charmness
   * \param b Bottomness
   * \param B Baryon number
   */
  QuantumNumbers(const FourVector& m,
                 const int q, const int i3,
                 const int s, const int c,
                 const int b, const int B)
        : momentum_(m),
          charge_(q),
          isospin3_(i3),
          strangeness_(s),
          charmness_(c),
          bottomness_(b),
          baryon_number_(B) {}

  /** construct QuantumNumbers collection from the conserved quantities
   * found in \p particles
   */
  QuantumNumbers(const Particles &particles) {
    count_conserved_values(particles);
  }

  /** sum up all conserved quantities from \p particles
   *
   * (sets everything to 0 at the beginning)
   */
  void count_conserved_values(const Particles &particles) {
    momentum_ = FourVector(0,0,0,0);
    charge_ = isospin3_ = strangeness_ = charmness_
            = bottomness_ = baryon_number_ = 0;
    for (const ParticleData &data : particles.data()) {
      momentum_      += data.momentum();
      charge_        += data.pdgcode().charge();
      isospin3_      += data.pdgcode().isospin3();
      strangeness_   += data.pdgcode().strangeness();
      charmness_     += data.pdgcode().charmness();
      bottomness_    += data.pdgcode().bottomness();
      baryon_number_ += data.pdgcode().baryon_number();
    }
  }

  /** returns the total momentum four-vector \f$P^\mu = \sum_{i \in
   * \mbox{particles}} (E_i, \vec p_i)\f$
   *
   * \see QuantumNumbers::momentum_
   * \see ParticleData::momentum()
   */
  FourVector momentum() const {
    return momentum_;
  }
  /** returns the total charge \f$Q = \sum_{i \in \mbox{particles}}
   * q_i\f$
   *
   * \see QuantumNumbers::charge_
   * \see PdgCode::charge()
   */
  int charge() const {
    return charge_;
  }
  /** returns twice the total isospin-3 component (twice) \f$I = \sum_{i
   * \in \mbox{particles}} 2{I_3}_i\f$
   *
   * \see QuantumNumbers::isospin3_
   * \see PdgCode::isospin3()
   */
  int isospin3() const {
    return isospin3_;
  }
  /** returns the total strangeness \f$S = \sum_{i \in \mbox{particles}}
   * S_i\f$
   *
   * \see QuantumNumbers::strangeness_
   * \see PdgCode::strangeness()
   */
  int strangeness() const {
    return strangeness_;
  }
  /** returns the total charm \f$C = \sum_{i \in \mbox{particles}}
   * C_i\f$
   *
   * \see QuantumNumbers::charmness_
   * \see PdgCode::charmness()
   */
  int charmness() const {
    return charmness_;
  }
  /** returns the total bottom \f$b = \sum_{i \in \mbox{particles}}
   * b_i\f$
   *
   * \see QuantumNumbers::bottomness_
   * \see PdgCode::bottomness()
   */
  int bottomness() const {
    return bottomness_;
  }
  /** returns the total baryon number \f$B = \sum_{i \in
   * \mbox{particles}} B_i\f$
   *
   * \see QuantumNumbers::baryon_number_
   * \see PdgCode::baryon_number()
   */
  int baryon_number() const {
    return baryon_number_;
  }

  /** returns true if all members compare true.
   *
   * In the comparison of FourVectors, a little leeway is built-in, so
   * this does not rely on an exact comparison of floating point values.
   *
   * \see FourVector::operator==
   *
   */
  bool operator==(const QuantumNumbers& rhs) const {
    return (momentum_      == rhs.momentum_
         && charge_        == rhs.charge_
         && isospin3_      == rhs.isospin3_
         && strangeness_   == rhs.strangeness_
         && charmness_     == rhs.charmness_
         && bottomness_    == rhs.bottomness_
         && baryon_number_ == rhs.baryon_number_);
  }
  /// logical complement of QuantumNumbers::operator==
  bool operator!=(const QuantumNumbers& rhs) const {
    return !(*this == rhs);
  }

  /** entry-wise subtract of two sets of QuantumNumbers
   *
   * If everything is conserved, all entries of the result should be
   * zero.
   *
   */
  QuantumNumbers operator-(const QuantumNumbers& rhs) const {
    return {momentum_      - rhs.momentum_,
            charge_        - rhs.charge_,
            isospin3_      - rhs.isospin3_,
            strangeness_   - rhs.strangeness_,
            charmness_     - rhs.charmness_,
            bottomness_    - rhs.bottomness_,
            baryon_number_ - rhs.baryon_number_};
  }

  /** checks if the current particle list has still the same values and
   * reports about differences.
   *
   * \see QuantumNumbers::report_deviations(const QuantumNumbers&) const
   */
  std::string report_deviations(const Particles& particles) const {
    QuantumNumbers current_values(particles);
    return report_deviations(current_values);
  }

  /** reports on deviations between two QuantumNumbers collections.
   *
   * \param rhs other QuantumNumbers collection
   * \return a string with information about the differences.
   *
   * If there are no differences, the returned string is empty; else, a
   * descriptive warning message is returned, e.g.
   *
   * \code
   * Conservation law violations detected (old vs. new)
   * Deviation in Charge:
   *  164 vs. 163
   * Deviation in Isospin 3:
   *  -88 vs -90
   * \endcode
   *
   */
  std::string report_deviations(const QuantumNumbers& rhs) const {
    if (rhs == *this) {
      return "";
    }
    std::stringstream error_msg;
    error_msg << "Conservation law violations detected (old vs. new)\n";
    if (momentum_ != rhs.momentum_) {
      error_msg << "Deviation in Four-Momentum:\n" << std::scientific;
    }
    // programmer's note: here, I'd like to simultaneously loop over an
    // integer (for the output; so that we know which component is
    // faulty) and both the current and rhs's momentum four-vector. If
    // there is a better way to do this, feel free to implement.
    //
    // I chose mu < 4 as the breaking condition out of the vague feeling
    // that comparing integers may be faster than accessing the
    // iterators.
    int mu = 0;
    for (auto here_iter = momentum_.cbegin(),
              rhs_iter = rhs.momentum_.cbegin();
         mu < 4;
         ++here_iter, ++rhs_iter, ++mu) {
      if (! almost_equal(*here_iter, *rhs_iter)) {
        error_msg << " P_" << mu << ": " << *here_iter << " vs. " << *rhs_iter
                  << "; Δ = " << (*here_iter - *rhs_iter) << "\n";
      }
     }
    if (charge_ != rhs.charge_) {
      error_msg << "Deviation in Charge:\n " << charge_ << " vs. "
                << rhs.charge_ << "\n";
    }
    if (isospin3_ != rhs.isospin3_) {
      error_msg << "Deviation in Isospin 3:\n " << isospin3_ << " vs ."
                << rhs.isospin3_ << "\n";
    }
    if (strangeness_ != rhs.strangeness_) {
      error_msg << "Deviation in Strangeness:\n " << strangeness_ << " vs. "
                << rhs.strangeness_ << "\n";
    }
    if (charmness_ != rhs.charmness_) {
      error_msg << "Deviation in Charmness:\n " << charmness_ << " vs. "
                << rhs.charmness_ << "\n";
    }
    if (bottomness_ != rhs.bottomness_) {
      error_msg << "Deviation in Bottomness:\n " << bottomness_ << " vs. "
                << rhs.bottomness_ << "\n";
    }
    if (baryon_number_ != rhs.baryon_number_) {
      error_msg << "Deviation in Baryon Number:\n " << baryon_number_ << " vs. "
                << rhs.baryon_number_ << "\n";
    }
    return error_msg.str();
  }

 private:
  /** total momentum four-vector
   *
   * \see QuantumNumbers::momentum()
   */
  FourVector momentum_;
  /** total charge
   *
   * \see QuantumNumbers::charge()
   */
  int charge_;
  /** total isospin-3
   *
   * \see QuantumNumbers::isospin3()
   */
  int isospin3_;
  /** total strangeness
   *
   * \see QuantumNumbers::strangeness()
   */
  int strangeness_;
  /** total charm
   *
   * \see QuantumNumbers::charmness()
   */
  int charmness_;
  /** total bottom
   *
   * \see QuantumNumbers::bottomness()
   */
  int bottomness_;
  /** total baryon number
   *
   * \see QuantumNumbers::baryon_number()
   */
  int baryon_number_;
};

}  // namespace SMASH

#endif  // SRC_INCLUDE_PDGCODE_H_
