/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PARTICLETYPE_H_
#define SRC_INCLUDE_PARTICLETYPE_H_

#include <string>
#include "include/pdgcode.h"

namespace Smash {

/**
 * Particle type contains the static properties of a particle
 *
 * SMASH reads in first the list of particles with their properties
 * and they don't change from this time on.
 * They are looked up according to their PDG code.
 */
class ParticleType {
 public:
  /// Use improbable values for default constructor
  ParticleType()
      : name_("unknown"),
        mass_(-1),
        width_(-1),
        pdgcode_(0x0),
        isospin_(0),
        charge_(0) {}
  /// Explicit constructor
  ParticleType(std::string n, float m, float w, PdgCode id)
      : name_(n),
        mass_(m),
        width_(w),
        pdgcode_(id),
        isospin_(pdgcode_.isospin_total()),
        charge_(pdgcode_.charge())
        {}
  /// set particle type
  inline void set(const std::string &n, float m, float w, PdgCode id);
  /// return particle name
  inline std::string name(void) const;
  /// return particle mass
  inline float mass(void) const;
  /// return particle width
  inline float width(void) const;
  /// return particle pdgcode
  inline PdgCode pdgcode(void) const;
  /// Isospin is 2 * particle data book value
  inline int isospin(void) const;
  /// return particle charge
  inline int charge(void) const;
  /// Spin is 2 * particle data book value
  inline int spin(void) const;

 private:
  /// name of the particle
  /// \todo This variable is only used for debug output. Maybe `ifdef` it out
  ///       for non-debug builds to save the memory?
  std::string name_;
  /// mass of the particle
  float mass_;
  /// width of the particle
  float width_;
  /// PDG Code of the particle
  PdgCode pdgcode_;
  /** twice the isospin of the particle
   *
   * This is filled automatically from pdgcode_.
   */
  int isospin_;
  /** charge of the particle
   *
   * This is filled automatically from pdgcode_.
   */
  int charge_;
};

inline void ParticleType::set(const std::string &NAME, float MASS,
     float WIDTH, PdgCode ID) {
  mass_ = MASS;
  width_ = WIDTH;
  pdgcode_ = PdgCode(ID);
  name_ = NAME;
  isospin_ = pdgcode_.isospin_total();
  charge_ = pdgcode_.charge();
}

inline int ParticleType::charge(void) const {
  return charge_;
}

inline int ParticleType::isospin(void) const {
  return isospin_;
}

inline float ParticleType::mass(void) const {
  return mass_;
}

inline std::string ParticleType::name(void) const {
  return name_;
}

inline PdgCode ParticleType::pdgcode(void) const {
  return pdgcode_;
}

inline int ParticleType::spin(void) const {
  return pdgcode_.spin();
}

inline float ParticleType::width(void) const {
  return width_;
}

}  // namespace Smash

#endif  // SRC_INCLUDE_PARTICLETYPE_H_
