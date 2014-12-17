/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PARTICLEDATA_H_
#define SRC_INCLUDE_PARTICLEDATA_H_

#include "forwarddeclarations.h"
#include "fourvector.h"
#include "particletype.h"
#include "pdgcode.h"

namespace Smash {

/**
 * \ingroup data
 *
 * ParticleData contains the dynamic information of a certain particle.
 *
 * Each particle has its momentum, position and other relevant physical
 * data entry.
 */
class ParticleData {
 public:
  /**
   * Create a new particle with the given \p particle_type and optionally a
   * specific \p unique_id.
   *
   * All other values are initialized to improbable values.
   */
  explicit ParticleData(const ParticleType &particle_type, int unique_id = -1)
      : id_(unique_id), type_(&particle_type) {}

  inline int id(void) const;
  inline void set_id(const int id);
  inline PdgCode pdgcode(void) const;

  // convenience accessors to PdgCode:
  /// \copydoc PdgCode::is_hadron
  bool is_hadron() const { return type_->is_hadron(); }

  /// \copydoc PdgCode::is_baryon
  bool is_baryon() const { return pdgcode().is_baryon(); }

  /** Returns the particle's pole mass ("on-shell"). */
  float pole_mass() const { return type_->mass(); }
  /** Returns the particle's effective mass
   * (as determined from the 4-momentum, possibly "off-shell"). */
  float effective_mass() const { return momentum().abs(); }

  /**
   * Return the ParticleType object associated to this particle.
   */
  const ParticleType &type() const { return *type_; }

  inline int id_process(void) const;
  inline void set_id_process(const int id);
  inline const FourVector &momentum(void) const;
  inline void set_4momentum(const FourVector &momentum_vector);
  inline void set_4momentum(double mass, const ThreeVector &mom);
  inline void set_4momentum(double mass, double px, double py, double pz);
  inline void set_3momentum(const ThreeVector &mom);
  inline const FourVector &position(void) const;
  inline void set_4position(const FourVector &pos);
  inline void set_3position(const ThreeVector &pos);
  /// get the velocity 3-vector
  inline ThreeVector velocity (void) const { return momentum_.threevec() / momentum_.x0(); }

  /**
   * Returns the inverse of the gamma factor from the current velocity of the
   * particle.
   *
   * \f[\frac{1}{\gamma}=\sqrt{1-v^2}\f]
   *
   * This functions is more efficient than calculating the gamma factor from
   * \ref velocity, since the \ref velocity function must execute three
   * divisions (for every space component of the momentum vector).
   */
  double inverse_gamma() const {
    return std::sqrt(1. -
                     (momentum_.x1() * momentum_.x1() +
                      momentum_.x2() * momentum_.x2() +
                      momentum_.x3() * momentum_.x3()) /
                         (momentum_.x0() * momentum_.x0()));
  }

  /// do a Lorentz-boost
  inline void boost (const ThreeVector &v);

  /* overloaded operators */
  inline bool operator==(const ParticleData &a) const;
  inline bool operator<(const ParticleData &a) const;
  inline bool operator==(const int id_a) const;
  inline bool operator<(const int id_a) const;

 private:
  /// Each particle has a unique identifier
  int id_ = -1;

  /**
   * A reference to the ParticleType object for this particle (this contains
   * all the static information).
   */
  ParticleTypePtr type_ = nullptr;

  /// counter of the last collision/decay
  int id_process_ = -1;
  /// momenta of the particle: x0, x1, x2, x3 as E, px, py, pz
  FourVector momentum_;
  /// position in space: x0, x1, x2, x3 as t, x, y, z
  FourVector position_;
};

/// look up the id of the particle
inline int ParticleData::id(void) const {
  return id_;
}

/// set id of the particle
inline void ParticleData::set_id(const int i) {
  id_ = i;
}

/// look up the pdgcode of the particle
inline PdgCode ParticleData::pdgcode(void) const {
  return type_->pdgcode();
}

/// look up the id of the collision process
inline int ParticleData::id_process(void) const {
  return id_process_;
}

/// set the id of the collision process
inline void ParticleData::set_id_process(const int process_id) {
  id_process_ = process_id;
}

/// return the particle 4-momentum
inline const FourVector &ParticleData::momentum(void) const {
  return momentum_;
}

/// set particle 4-momentum directly
inline void ParticleData::set_4momentum(const FourVector &momentum_vector) {
  momentum_ = momentum_vector;
}

/** set particle 4-momentum
 *
 * \param[in] new_mass the mass of the particle
 * \param[in] mom the three-momentum of the particle
 */
inline void ParticleData::set_4momentum(double new_mass, const ThreeVector &mom) {
  momentum_ = FourVector(std::sqrt(new_mass * new_mass + mom * mom), mom);
}

/** set particle 4-momentum
 *
 * \param[in] new_mass the mass of the particle
 * \param[in] px x-component of the momentum
 * \param[in] py y-component of the momentum
 * \param[in] pz z-component of the momentum
 */
inline void ParticleData::set_4momentum(double new_mass, double px, double py,
                                       double pz) {
  momentum_ = FourVector(
      std::sqrt(new_mass * new_mass + px * px + py * py + pz * pz), px, py, pz);
}

/// set particle 3-momentum
inline void ParticleData::set_3momentum(const ThreeVector &mom) {
  momentum_ = FourVector(momentum_.x0(), mom);
}

/// particle position in Minkowski space
inline const FourVector &ParticleData::position(void) const {
  return position_;
}

/// set the particle position directly
inline void ParticleData::set_4position(const FourVector &pos) {
  position_ = pos;
}

/// set the particle 3-position only (time component is not changed)
inline void ParticleData::set_3position(const ThreeVector &pos) {
  position_ = FourVector(position_.x0(), pos);
}

inline void ParticleData::boost (const ThreeVector &v)
{
  set_4momentum(momentum_.LorentzBoost(v));
  // we also need to boost the position
  set_4position(position_.LorentzBoost(v));
}

/// check if the particles are identical
inline bool ParticleData::operator==(const ParticleData &a) const {
  return this->id_ == a.id_;
}

/// sort particles along their id
inline bool ParticleData::operator<(const ParticleData &a) const {
  return this->id_ < a.id_;
}

/// check if the particles are identical to a given id
inline bool ParticleData::operator==(const int id_a) const {
  return this->id_ == id_a;
}

/// sort particles along to given id
inline bool ParticleData::operator<(const int id_a) const {
  return this->id_ < id_a;
}

/**\ingroup logging
 * Writes the state of the particle to the output stream.
 */
std::ostream &operator<<(std::ostream &s, const ParticleData &p);

/** \ingroup logging
 * Writes a compact overview over the particles in the \p particle_list argument
 * to the stream.
 */
std::ostream &operator<<(std::ostream &out, const ParticleList &particle_list);

}  // namespace Smash

#endif  // SRC_INCLUDE_PARTICLEDATA_H_
