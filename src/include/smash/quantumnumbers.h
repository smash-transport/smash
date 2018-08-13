/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_QUANTUMNUMBERS_H_
#define SRC_INCLUDE_QUANTUMNUMBERS_H_

#include <string>

#include "particles.h"

namespace smash {

/**
 * \ingroup data
 *
 * A container for storing conserved values.
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
  /// Construct QuantumNumbers collection with all fields 0.
  QuantumNumbers()
      : momentum_(0., 0., 0., 0.),
        charge_(0),
        isospin3_(0),
        strangeness_(0),
        charmness_(0),
        bottomness_(0),
        baryon_number_(0) {}

  /**
   * \return Constructed QuantumNumbers.
   * \param[in] m Momentum FourVector [GeV]
   * \param[in] q Charge
   * \param[in] i3 Isospin
   * \param[in] s Strangeness
   * \param[in] c Charmness
   * \param[in] b Bottomness
   * \param[in] B Baryon number
   */
  QuantumNumbers(const FourVector& m, const int q, const int i3, const int s,
                 const int c, const int b, const int B)
      : momentum_(m),
        charge_(q),
        isospin3_(i3),
        strangeness_(s),
        charmness_(c),
        bottomness_(b),
        baryon_number_(B) {}

  /**
   * Construct QuantumNumbers collection from the conserved quantities
   * found in a set of particles.
   * \param[in] particles set of particles for which
   *            quantum numbers are calculated and constructed
   * \return Constructed object.
   */
  explicit QuantumNumbers(const Particles& particles) : QuantumNumbers() {
    for (const ParticleData& data : particles) {
      add_values(data);
    }
  }

  /**
   * Construct QuantumNumbers collection from a particle list.
   * \param[in] part list of particles for which
   *            quantum numbers are calculated and constructed
   * \return Constructed object.
   */
  explicit QuantumNumbers(const ParticleList& part) : QuantumNumbers() {
    for (const auto& p : part) {
      add_values(p);
    }
  }

  /**
   * Add the quantum numbers of a single particle to the collection.
   * \param[in] p particle whose quantum number is added to the collection
   */
  void add_values(const ParticleData& p) {
    momentum_ += p.momentum();
    charge_ += p.pdgcode().charge();
    isospin3_ += p.pdgcode().isospin3();
    strangeness_ += p.pdgcode().strangeness();
    charmness_ += p.pdgcode().charmness();
    bottomness_ += p.pdgcode().bottomness();
    baryon_number_ += p.pdgcode().baryon_number();
  }

  /**
   * \return The total momentum four-vector.
   * \f$P^\mu = \sum_{i \in \mbox{particles}} (E_i, \vec p_i)\f$ [GeV]
   *
   * \see QuantumNumbers::momentum_
   * \see ParticleData::momentum()
   */
  FourVector momentum() const { return momentum_; }

  /**
   * \return The total electric charge.
   * \f$Q = \sum_{i \in \mbox{particles}} q_i\f$
   *
   * \see QuantumNumbers::charge_
   * \see PdgCode::charge()
   */
  int charge() const { return charge_; }
  /**
   * \return Twice the total isospin-3 component.
   * \f$I = \sum_{i \in \mbox{particles}} 2{I_3}_i\f$
   *
   * \see QuantumNumbers::isospin3_
   * \see PdgCode::isospin3()
   */
  int isospin3() const { return isospin3_; }
  /**
   * \return The total strangeness.
   * \f$S = \sum_{i \in \mbox{particles}} S_i\f$
   *
   * \see QuantumNumbers::strangeness_
   * \see PdgCode::strangeness()
   */
  int strangeness() const { return strangeness_; }
  /**
   * \return The total charm.
   * \f$C = \sum_{i \in \mbox{particles}} C_i\f$
   *
   * \see QuantumNumbers::charmness_
   * \see PdgCode::charmness()
   */
  int charmness() const { return charmness_; }
  /**
   * \return The total bottom.
   * \f$b = \sum_{i \in \mbox{particles}} b_i\f$
   *
   * \see QuantumNumbers::bottomness_
   * \see PdgCode::bottomness()
   */
  int bottomness() const { return bottomness_; }
  /**
   * \return The total baryon number.
   * \f$B = \sum_{i \in \mbox{particles}} B_i\f$
   *
   * \see QuantumNumbers::baryon_number_
   * \see PdgCode::baryon_number()
   */
  int baryon_number() const { return baryon_number_; }

  /**
   * \return true if all members compare true.
   * \param rhs Right-hand side.
   *
   * In the comparison of FourVectors, a little leeway is built-in, so
   * this does not rely on an exact comparison of floating point values.
   *
   * \see FourVector::operator==
   */
  bool operator==(const QuantumNumbers& rhs) const {
    // clang-format off
    return (momentum_ == rhs.momentum_ &&
            charge_ == rhs.charge_ &&
            isospin3_ == rhs.isospin3_ &&
            strangeness_ == rhs.strangeness_ &&
            charmness_ == rhs.charmness_ &&
            bottomness_ == rhs.bottomness_ &&
            baryon_number_ == rhs.baryon_number_);
    // clang-format on
  }
  /// Logical complement of QuantumNumbers::operator==.
  bool operator!=(const QuantumNumbers& rhs) const { return !(*this == rhs); }

  /**
   * \return Entry-wise difference of two sets of QuantumNumbers.
   * \param rhs Right-hand side.
   *
   * If everything is conserved, all entries of the result should be
   * zero.
   */
  QuantumNumbers operator-(const QuantumNumbers& rhs) const {
    return {momentum_ - rhs.momentum_,          charge_ - rhs.charge_,
            isospin3_ - rhs.isospin3_,          strangeness_ - rhs.strangeness_,
            charmness_ - rhs.charmness_,        bottomness_ - rhs.bottomness_,
            baryon_number_ - rhs.baryon_number_};
  }

  /**
   * Checks if the current particle list has still the same values and
   * reports about differences.
   * \param[in] particles Set of particles whose quantum number is compared
   * \return String reporting the deviations.
   *
   * \see QuantumNumbers::report_deviations(const QuantumNumbers&) const
   */
  std::string report_deviations(const Particles& particles) const {
    QuantumNumbers current_values(particles);
    return report_deviations(current_values);
  }

  /**
   * Reports on deviations between two QuantumNumbers collections.
   *
   * \param[in] rhs Other QuantumNumbers collection.
   * \return A string with information about the differences.
   *
   * If there are no differences, the returned string is empty; else, a
   * descriptive warning message is returned, e.g.
   *
   * \code
   * Conservation law violations detected (old vs. new)
   * Deviation in Charge:
   *  164 vs. 163
   * Deviation in Isospin 3:
   *  -88 vs. -90
   * \endcode
   *
   */
  std::string report_deviations(const QuantumNumbers& rhs) const;

 private:
  /**
   * Total momentum four-vector [GeV].
   *
   * \see QuantumNumbers::momentum()
   */
  FourVector momentum_;

  /**
   * Total charge.
   *
   * \see QuantumNumbers::charge()
   */
  int charge_;
  /**
   * Total isospin-3.
   *
   * \see QuantumNumbers::isospin3()
   */
  int isospin3_;
  /**
   * Total strangeness.
   *
   * \see QuantumNumbers::strangeness()
   */
  int strangeness_;
  /**
   * Total charm.
   *
   * \see QuantumNumbers::charmness()
   */
  int charmness_;
  /**
   * Total bottom.
   *
   * \see QuantumNumbers::bottomness()
   */
  int bottomness_;
  /**
   * Total baryon number.
   *
   * \see QuantumNumbers::baryon_number()
   */
  int baryon_number_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_QUANTUMNUMBERS_H_
