/*
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_DECAYTYPE_H_
#define SRC_INCLUDE_DECAYTYPE_H_

#include <memory>
#include <vector>

#include "forwarddeclarations.h"
#include "particletype.h"
#include "tabulation.h"

namespace smash {

/**
 * DecayType is the abstract base class for all decay types.
 */
class DecayType {
 public:
  /**
   * Construct a \ref DecayType.
   *
   * \param[in] part_types Final-state particles of the decay.
   * \param[in] l Angular momentum of the decay.
   * \return The constructed object.
   */
  DecayType(ParticleTypePtrList part_types, int l)
      : particle_types_(part_types), L_(l) {}
  /**
   * Virtual Destructor.
   *
   * The declaration of the destructor is necessary to make it virtual.
   */
  virtual ~DecayType() = default;
  /// \return the number of particles in the final state
  virtual unsigned int particle_number() const = 0;
  /**
   * \return if the final state consists of the given particle list.
   *
   * \param[in] list Final state particle types to be checked.
   */
  virtual bool has_particles(ParticleTypePtrList list) const = 0;
  /**
   * \return if this decay type has the right mother
   * (most decays do not depend on the mother type).
   *
   * \param[in] mother Particle type to be checked.
   */
  virtual bool has_mother(ParticleTypePtr mother) const {
    SMASH_UNUSED(mother);
    return true;
  }
  /// \return the particle types associated with this branch.
  const ParticleTypePtrList &particle_types() const { return particle_types_; }
  /// \return the angular momentum of this branch.
  inline int angular_momentum() const { return L_; }
  /**
   * \return the mass-dependent width of the decay.
   *
   * \param[in] m0 Pole mass of the decaying particle [GeV].
   * \param[in] G0 Partial width at the pole mass [GeV].
   * \param[in] m Actual mass of the decaying particle [GeV].
   */
  virtual double width(double m0, double G0, double m) const = 0;
  /**
   * \return The mass-dependent in-width for a resonance formation process.
   *
   * \param[in] m0 Pole mass of the produced resonance [GeV].
   * \param[in] G0 Partial width at the pole mass [GeV].
   * \param[in] m Actual mass of the produced resonance [GeV].
   * \param[in] m1 Actual mass of the first incoming particle [GeV].
   * \param[in] m2 Actual mass of the second incoming particle [GeV].
   */
  virtual double in_width(double m0, double G0, double m, double m1,
                          double m2) const = 0;

 protected:
  /// final-state particles of the decay
  ParticleTypePtrList particle_types_;
  /// angular momentum of the decay
  int L_;
};

/**
 * TwoBodyDecay represents a decay type with two final-state particles.
 */
class TwoBodyDecay : public DecayType {
 public:
  /**
   * Construct a \ref TwoBodyDecay.
   *
   * \param[in] part_types Final-state particles of the decay.
   * \param[in] l Angular momentum of the decay.
   * \return The constructed object.
   */
  TwoBodyDecay(ParticleTypePtrList part_types, int l);
  unsigned int particle_number() const override;
  bool has_particles(ParticleTypePtrList list) const override;
  /// \return The kinematic energy threshold of the decay in GeV.
  double threshold() const {
    return particle_types_[0]->min_mass_kinematic() +
           particle_types_[1]->min_mass_kinematic();
  }

 protected:
  /**
   * This is a virtual helper method which is used to write the width as
   * Gamma(m) = Gamma_0 * rho(m) / rho(m_0). This ensures that the width is
   * properly normalized at the pole mass to Gamma(m_0) = Gamma_0.
   * By default rho simply equals one, which corresponds to a constant width.
   *
   * \param[in] mass Resonance mass of the decay.
   * \return Width of the decay at the given \p mass.
   * */
  virtual double rho(double mass) const {
    SMASH_UNUSED(mass);
    return 1.;
  }
};

/**
 * TwoBodyDecayStable represents a decay type with two stable final-state
 * particles.
 */
class TwoBodyDecayStable : public TwoBodyDecay {
 public:
  /**
   * Construct a \ref TwoBodyDecayStable.
   *
   * \param[in] part_types Final-state particles of the decay.
   * \param[in] l Angular momentum of the decay.
   * \return The constructed object.
   */
  TwoBodyDecayStable(ParticleTypePtrList part_types, int l);

  double width(double m0, double G0, double m) const override;

  double in_width(double m0, double G0, double m, double m1,
                  double m2) const override;

 protected:
  double rho(double m) const override;
};

/**
 * TwoBodyDecaySemistable represents a decay type with two final-state
 * particles, one of which is stable and the other is unstable.
 */
class TwoBodyDecaySemistable : public TwoBodyDecay {
 public:
  /**
   * Construct a \ref TwoBodyDecaySemistable.
   *
   * \param[in] part_types Final-state particles of the decay.
   * \param[in] l Angular momentum of the decay.
   * \return The constructed object.
   */
  TwoBodyDecaySemistable(ParticleTypePtrList part_types, int l);

  double width(double m0, double G0, double m) const override;

  double in_width(double m0, double G0, double m, double m1,
                  double m2) const override;

 protected:
  double rho(double m) const override;

  /**
   * \return the cutoff parameter Λ for semi-stable decays,
   * given the types of the daughter particles.
   *
   * For the values used in GiBUU, see \iref{Buss:2011mx}, eq. (175).
   * For the original values used by M. Post, see table 1 in \iref{Post:2003hu}.
   *
   * We mostly stick to the GiBUU values, but use a different value for the ρπ
   * decay, in order to avoid secondary bumps in the ω spectral function and
   * achieve a better normalization. In contrast to smash, GiBUU does not have
   * an ω → ρ π decay.
   */
  double get_Lambda();

  /// Cut-off parameter Λ for semi-stable decays.
  double Lambda_;

  /// Tabulation of the resonance integrals.
  mutable std::unique_ptr<Tabulation> tabulation_;
};

/**
 * TwoBodyDecayUnstable represents a decay type with two unstable final-state
 * particles.
 */
class TwoBodyDecayUnstable : public TwoBodyDecay {
 public:
  /**
   * Construct a \ref TwoBodyDecayUnstable.
   *
   * \param[in] part_types Final-state particles of the decay.
   * \param[in] l Angular momentum of the decay.
   * \return The constructed object.
   */
  TwoBodyDecayUnstable(ParticleTypePtrList part_types, int l);
  double width(double m0, double G0, double m) const override;
  double in_width(double m0, double G0, double m, double m1,
                  double m2) const override;

 protected:
  double rho(double m) const override;
  /**
   * \return the cut-off parameter Λ for unstable decays,
   * given the types of the daughter particles.
   */
  double get_Lambda();

  /// Cut-off parameter Λ for unstable decays.
  double Lambda_;

  /// Tabulation of the resonance integrals.
  mutable std::unique_ptr<Tabulation> tabulation_;
};

/**
 * TwoBodyDecayDilepton represents a decay with a lepton and its antilepton
 * as the final-state particles.
 */
class TwoBodyDecayDilepton : public TwoBodyDecayStable {
 public:
  /**
   * Construct a \ref TwoBodyDecayDilepton.
   *
   * \param[in] part_types Final-state particles of the decay.
   * \param[in] l Angular momentum of the decay.
   * \return The constructed object.
   */
  TwoBodyDecayDilepton(ParticleTypePtrList part_types, int l);

  double width(double m0, double G0, double m) const override;
};

/**
 * ThreeBodyDecay represents a decay type with three final-state particles.
 */
class ThreeBodyDecay : public DecayType {
 public:
  /**
   * Construct a \ref ThreeBodyDecay.
   *
   * \param[in] part_types Final-state particles of the decay.
   * \param[in] l Angular momentum of the decay.
   * \return The constructed object.
   */
  ThreeBodyDecay(ParticleTypePtrList part_types, int l);

  unsigned int particle_number() const override;
  bool has_particles(ParticleTypePtrList list) const override;
  double width(double m0, double G0, double m) const override;
  double in_width(double m0, double G0, double m, double m1,
                  double m2) const override;
};

/**
 * ThreeBodyDecayDilepton represents a decay type with three final-state
 * particles, two of which are leptons.
 */
class ThreeBodyDecayDilepton : public ThreeBodyDecay {
 public:
  /**
   * Construct a \ref ThreeBodyDecayDilepton.
   *
   * \param[in] mother Type of the mother particle.
   * \param[in] part_types Final-state particles of the decay.
   * \param[in] l Angular momentum of the decay.
   * \return The constructed object.
   */
  ThreeBodyDecayDilepton(ParticleTypePtr mother, ParticleTypePtrList part_types,
                         int l);

  bool has_mother(ParticleTypePtr mother) const override;

  /**
   * Get the mass-differential width \f$ d\Gamma / dm \f$ for a dilepton Dalitz
   * decay, where \f$ m \f$ is the invariant mass of the lepton pair.
   *
   * This differential width is used directly for the dilepton shining weights.
   * It is calculated according to \iref{Weil:2013mya}, eq. (30)-(36).
   *
   * \todo{Cite SMASH dilepton paper.}
   *
   * \param[in] m_par Mass of the parent.
   * \param[in] m_l Mass of the lepton species.
   * \param[in] m_dil Invariant mass of the dilepton pair [GeV].
   * \param[in] m_other Mass of the third, non-leptonic, particle.
   * \param[in] other Type of the third  particle.
   * \param[in] t Type of the parent particle.
   * \return Mass differential width.
   */
  static double diff_width(double m_par, double m_l, double m_dil,
                           double m_other, ParticleTypePtr other,
                           ParticleTypePtr t);
  double width(double m0, double G0, double m) const override;

 protected:
  /// Tabulation of the resonance integrals.
  mutable std::unique_ptr<Tabulation> tabulation_;

  /// Type of the mother particle.
  ParticleTypePtr mother_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_DECAYTYPE_H_
