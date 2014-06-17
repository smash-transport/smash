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

  /// \copydoc ParticleType::mass
  float mass() const { return type_->mass(); }

  /**
   * Return the ParticleType object associated to this particle.
   *
   * \todo Remove the need for the Particles parameter.
   */
  const ParticleType &type() const { return *type_; }

  inline int id_process(void) const;
  inline void set_id_process(const int id);
  inline double collision_time(void) const;
  inline void set_collision_time(const double &collision_time);
  inline void set_collision(const double &collision_time);
  inline void set_collision_past(const int process_id);
  inline const FourVector &momentum(void) const;
  inline void set_momentum(const FourVector &momentum_vector);
  inline void set_momentum(double mass, const ThreeVector &mom);
  inline void set_momentum(double mass, double px, double py, double pz);
  inline const FourVector &position(void) const;
  inline void set_position(const FourVector &position);
  /// get the velocity 3-vector
  inline ThreeVector velocity (void) const { return momentum_.threevec() / momentum_.x0(); }
  /// do a Lorentz-boost
  inline void boost (const ThreeVector &v);
  /// get the full decay width (mass-dependent!) of a particular particle
  float total_width() const;

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
  const ParticleType *type_ = nullptr;

  /// counter of the last collision/decay
  int id_process_ = -1;
  /// collision time
  double collision_time_ = 0.;
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

/// return collision time
inline double ParticleData::collision_time(void) const {
  return collision_time_;
}

/// set possible collision time
inline void ParticleData::set_collision_time(const double &collision_t) {
  collision_time_ = collision_t;
}

/// set possible collision data
inline void ParticleData::set_collision(const double &collision_t) {
  collision_time_ = collision_t;
}

/// set happened collision data
inline void ParticleData::set_collision_past(const int id_counter) {
  collision_time_ = 0.0;
  id_process_ = id_counter;
}

/// return the particle four momentum
inline const FourVector &ParticleData::momentum(void) const {
  return momentum_;
}

/// set particle four momentum directly
inline void ParticleData::set_momentum(const FourVector &momentum_vector) {
  momentum_ = momentum_vector;
}

/** set particle four momentum
 *
 * \param[in] new_mass the mass of the particle
 * \param[in] mom the three-momentum of the particle
 */
inline void ParticleData::set_momentum(double new_mass, const ThreeVector &mom) {
  momentum_ = FourVector(std::sqrt(new_mass * new_mass + mom * mom), mom);
}

/** set particle four momentum
 *
 * \param[in] new_mass the mass of the particle
 * \param[in] px x-component of the momentum
 * \param[in] py y-component of the momentum
 * \param[in] pz z-component of the momentum
 */
inline void ParticleData::set_momentum(double new_mass, double px, double py,
                                       double pz) {
  momentum_ = FourVector(
      std::sqrt(new_mass * new_mass + px * px + py * py + pz * pz), px, py, pz);
}

/// particle position in Minkowski space
inline const FourVector &ParticleData::position(void) const {
  return position_;
}

/// set the particle position directly
inline void ParticleData::set_position(const FourVector &pos) {
  position_ = pos;
}

inline void ParticleData::boost (const ThreeVector &v)
{
  set_momentum(momentum_.LorentzBoost(v));
  // TODO: do we actually need to boost the position?
  set_position(position_.LorentzBoost(v));
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

}  // namespace Smash

#endif  // SRC_INCLUDE_PARTICLEDATA_H_
