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
  virtual int particle_number() const = 0;
  // Check if the final state consists of the given two particles.
  virtual bool has_particles(const ParticleType &t_a,
                             const ParticleType &t_b) const = 0;
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
  int particle_number() const override;
  bool has_particles(const ParticleType &t_a,
                     const ParticleType &t_b) const override;

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
   * see also PhD thesis Effenberger, eq. (2.77).
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
  int particle_number() const override;
  bool has_particles(const ParticleType &, const ParticleType &) const override;
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
  ThreeBodyDecayDilepton(ParticleTypePtrList part_types, int l);
  /**
   * Get the differential width for dilepton dalitz decay. Because we use the
   * shining method, we do not need a partial width and can use the differential
   * width directly for the shinin weights. The differential width is calculated
   * according to PhD thesis Weil, eq. (30) - (36).
   */
  static float diff_width(float m_parent, float m_dil,
                   float m_other, PdgCode pdg);
  float width(float m0, float G0, float m) const override;
 protected:
  std::unique_ptr<Tabulation> tabulation_;
};


}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYTYPE_H_
