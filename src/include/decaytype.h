/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_DECAYTYPE_H_
#define SRC_INCLUDE_DECAYTYPE_H_

#include <vector>
#include <memory>

#include "forwarddeclarations.h"
#include "particletype.h"
#include "tabulation.h"

namespace Smash {

/**
 * DecayType is the abstract base class for all decay types.
 */
class DecayType {
 public:
  DecayType(ParticleTypePtrList part_types, int l)
           : particle_types_(part_types), L_(l) {}
  /** Virtual Destructor.
   * The declaration of the destructor is necessary to make it virtual.
   */
  virtual ~DecayType() = default;
  // Get the number of particles in the final state
  virtual unsigned int particle_number() const = 0;
  // Check if the final state consists of the given particle list
  virtual bool has_particles(ParticleTypePtrList list) const = 0;
  /* Check if this decay type has the right mother
   * (most decays do not depend on the mother type). */
  virtual bool has_mother(ParticleTypePtr) const { return true; }
  /// Return the particle types associated with this branch.
  const ParticleTypePtrList &particle_types() const {
    return particle_types_;
  }
  /// Get the angular momentum of this branch.
  inline int angular_momentum() const {
    return L_;
  }
  /**
   * Get the mass-dependent width of the decay.
   *
   * \param m0 Pole mass of the decaying particle [GeV].
   * \param G0 Partial width at the pole mass [GeV].
   * \param m Actual mass of the decaying particle [GeV].
   */
  virtual float width(float m0, float G0, float m) const = 0;
  /**
   * Get the mass-dependent in-width for a resonance formation process.
   *
   * \param m0 Pole mass of the produced resonance [GeV].
   * \param G0 Partial width at the pole mass [GeV].
   * \param m Actual mass of the produced resonance [GeV].
   * \param m1 Actual mass of the first incoming particle [GeV].
   * \param m2 Actual mass of the second incoming particle [GeV].
   */
  virtual float in_width(float m0, float G0, float m,
                         float m1, float m2) const = 0;


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
  TwoBodyDecay(ParticleTypePtrList part_types, int l);
  unsigned int particle_number() const override;
  bool has_particles(ParticleTypePtrList list) const override;

 protected:
  /* This is a virtual helper method which is used to write the width as
   * Gamma(m) = Gamma_0 * rho(m) / rho(m_0). This ensures that the width is
   * properly normalized at the pole mass to Gamma(m_0) = Gamma_0.
   * By default rho simply equals one, which corresponds to a constant width. */
  virtual float rho(float) const { return 1.; }
};


/**
 * TwoBodyDecayStable represents a decay type with two stable final-state particles.
 */
class TwoBodyDecayStable : public TwoBodyDecay {
 public:
  TwoBodyDecayStable(ParticleTypePtrList part_types, int l);
  /**
   * Get the mass-dependent width of a two-body decay into stable particles
   * according to \iref{Manley:1992yb}.
   *
   * \param m0 Pole mass of the decaying particle [GeV].
   * \param G0 Partial width at the pole mass [GeV].
   * \param m Actual mass of the decaying particle [GeV].
   */
  float width(float m0, float G0, float m) const override;
  float in_width(float m0, float G0, float m,
                 float m1, float m2) const override;
 protected:
  float rho(float m) const override;
};


/**
 * TwoBodyDecaySemistable represents a decay type with two final-state particles,
 * one of which is stable and the other is unstable.
 */
class TwoBodyDecaySemistable : public TwoBodyDecay {
 public:
  TwoBodyDecaySemistable(ParticleTypePtrList part_types, int l);
  /**
   * Get the mass-dependent width of a two-body decay into one stable and one
   * unstable particle according to \iref{Manley:1992yb}.
   *
   * \param m0 Pole mass of the decaying particle [GeV].
   * \param G0 Partial width at the pole mass [GeV].
   * \param m Actual mass of the decaying particle [GeV].
   */
  float width(float m0, float G0, float m) const override;
  /**
   * Get the mass-dependent in-width for a resonance formation process from one
   * stable and one unstable particle according to \iref{Manley:1992yb},
   * see also \iref{Effenberger:1999wlg}, eq. (2.77).
   *
   * \param m0 Pole mass of the produced resonance [GeV].
   * \param G0 Partial width at the pole mass [GeV].
   * \param m Actual mass of the produced resonance [GeV].
   * \param m1 Actual mass of the first incoming particle [GeV].
   * \param m2 Actual mass of the second incoming particle [GeV].
   */
  float in_width(float m0, float G0, float m,
                 float m1, float m2) const override;

 protected:
  float rho(float m) const override;
  /**
   * Determine the cutoff parameter Λ for semi-stable decays,
   * given the types of the daughter particles.
   *
   * For the values used in GiBUU, see \iref{Buss:2011mx}, eq. (175).
   * For the original values used by M. Post, see table 1 in \iref{Post:2003hu}.
   *
   * We mostly stick to the GiBUU values, but use a different value for the ρπ
   * decay, in order to avoid secondary bumps in the ω spectral function and
   * achieve a better normalization. In contrast to Smash, GiBUU does not have
   * an ω → ρ π decay.
   */
  float get_Lambda();
  float Lambda_;
  std::unique_ptr<Tabulation> tabulation_;
};


/**
 * TwoBodyDecayUnstable represents a decay type with two unstable final-state particles.
 */
class TwoBodyDecayUnstable : public TwoBodyDecay {
 public:
  TwoBodyDecayUnstable(ParticleTypePtrList part_types, int l);
  float width(float m0, float G0, float m) const override;
  float in_width(float m0, float G0, float m,
                 float m1, float m2) const override;
 protected:
  float rho(float m) const override;
  /**
   * Determine the cutoff parameter Λ for unstable decays,
   * given the types of the daughter particles.
   */
  float get_Lambda();
  float Lambda_;
  std::unique_ptr<Tabulation> tabulation_;
};

/**
 * TwoBodyDecayDilepton represents a decay with a lepton and it's antilepton
 * as the final state particles.
 */
class TwoBodyDecayDilepton : public TwoBodyDecayStable {
 public:
  TwoBodyDecayDilepton(ParticleTypePtrList part_types, int l);
  float width(float m0, float G0, float m) const override;
};

/**
 * ThreeBodyDecay represents a decay type with three final-state particles.
 */
class ThreeBodyDecay : public DecayType {
 public:
  ThreeBodyDecay(ParticleTypePtrList part_types, int l);
  unsigned int particle_number() const override;
  bool has_particles(ParticleTypePtrList list) const override;
  float width(float m0, float G0, float m) const override;
  float in_width(float m0, float G0, float m,
                 float m1, float m2) const override;
};

/**
 * ThreeBodyDecayDilepton represents a decay type with three final-state
 * particles. Two of them are a dilepton.
 */
class ThreeBodyDecayDilepton : public ThreeBodyDecay {
 public:
  ThreeBodyDecayDilepton(ParticleTypePtr mother,
                         ParticleTypePtrList part_types, int l);
  bool has_mother(ParticleTypePtr mother) const override;
  /**
   * Get the differential width for dilepton dalitz decay. Because we use the
   * shining method, we do not need a partial width and can use the differential
   * width directly for the shining weights. The differential width is
   * calculated according to \iref{Weil:2013mya}, eq. (30)-(36).
   */
  static float diff_width(float m_parent, float m_dil,
                          float m_other, PdgCode pdg);
  float width(float m0, float G0, float m) const override;
 protected:
  std::unique_ptr<Tabulation> tabulation_;
  ParticleTypePtr mother_;
};


}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYTYPE_H_
