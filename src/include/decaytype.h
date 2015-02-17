/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_DECAYTYPES_H_
#define SRC_INCLUDE_DECAYTYPES_H_ 

#include "forwarddeclarations.h"
#include "particletype.h"

namespace Smash {

/**
 * DecayType is the abstract base class for all decay types.
 */
class DecayType {
 public:
  DecayType(ParticleTypePtrList part_types, int l) : particle_types_(part_types),
                                                      L_(l) {}
  /** Virtual Destructor.
   * The declaration of the destructor is necessary to make it virtual.
   */
  virtual ~DecayType() = default;
  // Get the number of particles in the final state
  virtual int particle_number() const = 0;
  // Check if the final state consists of the given two particles.
  virtual bool has_particles (const ParticleType &t_a, const ParticleType &t_b) const = 0;
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
  virtual float width (float m0, float G0, float m) const = 0;
  /**
   * Get the mass-dependent in-width for a resonance formation process.
   *
   * \param m0 Pole mass of the produced resonance [GeV].
   * \param G0 Partial width at the pole mass [GeV].
   * \param m Actual mass of the produced resonance [GeV].
   * \param m1 Actual mass of the first incoming particle [GeV].
   * \param m2 Actual mass of the second incoming particle [GeV].
   */
  virtual float in_width (float m0, float G0, float m, float m1, float m2) const = 0;

 protected:
  // Evaluate rho_ab according to equ. (2.76) in Effenberger's thesis.
  virtual float rho (float m) const = 0;
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
  TwoBodyDecay (ParticleTypePtrList part_types, int l);
  int particle_number() const override;
  bool has_particles (const ParticleType &t_a, const ParticleType &t_b) const override;
};


/**
 * TwoBodyDecayStable represents a decay type with two stable final-state particles.
 */
class TwoBodyDecayStable : public TwoBodyDecay {
 public:
  TwoBodyDecayStable (ParticleTypePtrList part_types, int l);
  /**
   * Get the mass-dependent width of a two-body decay into stable particles
   * according to Manley/Saleski, Phys. Rev. D 45 (1992) 4002.
   *
   * \param m0 Pole mass of the decaying particle [GeV].
   * \param G0 Partial width at the pole mass [GeV].
   * \param m Actual mass of the decaying particle [GeV].
   */
  float width (float m0, float G0, float m) const override;
  float in_width (float m0, float G0, float m, float m1, float m2) const override;
 protected:
  float rho (float m) const override;
};


/**
 * TwoBodyDecaySemistable represents a decay type with two final-state particles,
 * one of which is stable and the other is unstable.
 */
class TwoBodyDecaySemistable : public TwoBodyDecay {
 public:
  TwoBodyDecaySemistable (ParticleTypePtrList part_types, int l);
  /**
   * Get the mass-dependent width of a two-body decay into one stable and one
   * unstable particle according to Manley/Saleski, Phys. Rev. D 45 (1992) 4002.
   *
   * \param m0 Pole mass of the decaying particle [GeV].
   * \param G0 Partial width at the pole mass [GeV].
   * \param m Actual mass of the decaying particle [GeV].
   * \param m_uns Actual mass of the unstable incoming particle [GeV].
   */
  float width (float m0, float G0, float m) const override;
  /**
   * Get the mass-dependent in-width for a resonance formation process from one
   * stable and one unstable particle according to Manley/Saleski (PRD45),
   * see also PhD thesis Effenberger, eq. (2.77).
   *
   * \param m0 Pole mass of the produced resonance [GeV].
   * \param G0 Partial width at the pole mass [GeV].
   * \param m Actual mass of the produced resonance [GeV].
   * \param m1 Actual mass of the first incoming particle [GeV].
   * \param m2 Actual mass of the second incoming particle [GeV].
   */
  float in_width (float m0, float G0, float m, float m1, float m2) const override;
 protected:
  float rho (float m) const override;
  float calc_rho (float m) const;
  std::vector<float> tabulation_;  // vector for storing tabulated values
  float M_min_, dM_;               // minimum mass and step size for tabulation
  float Lambda_;
  void init_tabulation(float range, unsigned int N);
};


/**
 * TwoBodyDecayUnstable represents a decay type with two unstable final-state particles.
 */
class TwoBodyDecayUnstable : public TwoBodyDecay {
 public:
  TwoBodyDecayUnstable (ParticleTypePtrList part_types, int l);
  float width (float m0, float G0, float m) const override;
  float in_width (float m0, float G0, float m, float m1, float m2) const override;
 protected:
  float rho (float m) const override;
};


/**
 * ThreeBodyDecay represents a decay type with three final-state particles.
 */
class ThreeBodyDecay : public DecayType {
 public:
  ThreeBodyDecay (ParticleTypePtrList part_types, int l);
  int particle_number() const override;
  bool has_particles (const ParticleType &, const ParticleType &) const override;
  float width (float m0, float G0, float m) const override;
  float in_width (float m0, float G0, float m, float m1, float m2) const override;
 protected:
  float rho (float m) const override;
};


}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYTYPES_H_
