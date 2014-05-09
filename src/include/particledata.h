/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PARTICLEDATA_H_
#define SRC_INCLUDE_PARTICLEDATA_H_

#include <math.h>
#include <vector>

#include "include/fourvector.h"
#include "include/pdgcode.h"

namespace Smash {

/**
 * ParticleData contains the dynamic information of a certain particle.
 *
 * Each particle has its momentum, position and other relevant physical
 * data entry.
 */
class ParticleData {
 public:
  /// Use improbable values for default constructor
  ParticleData() :id_(-1), pdgcode_(0x0), id_process_(-1),
    collision_time_(0.0) {}
  /// Use improbable values for constructor
  explicit ParticleData(int i) :id_(i), pdgcode_(0x0),
    id_process_(-1), collision_time_(0.0) {}
  inline int id(void) const;
  inline void set_id(int id);
  inline PdgCode pdgcode(void) const;
  inline void set_pdgcode(PdgCode pdgcode);
  inline int id_process(void) const;
  inline void set_id_process(int id);
  inline double collision_time(void) const;
  inline void set_collision_time(const double &collision_time);
  inline void set_collision(const double &collision_time);
  inline void set_collision_past(int process_id);
  inline const FourVector &momentum(void) const;
  inline void set_momentum(const FourVector &momentum_vector);
  inline void set_momentum(const double &mass, const double &px,
                           const double &py, const double &pz);
  inline const FourVector &position(void) const;
  inline void set_position(const FourVector &position);
  inline void set_position(const double &x0, const double &x1,
                           const double &x2, const double &x3);
  /// get the x veolcity
  inline double velocity_x(void) { return momentum().x1() / momentum().x0(); }
  /// get the y veolcity
  inline double velocity_y(void) { return momentum().x2() / momentum().x0(); }
  /// get the z veolcity
  inline double velocity_z(void) { return momentum().x3() / momentum().x0(); }
  /* overloaded operators */
  inline bool operator==(const ParticleData &a) const;
  inline bool operator<(const ParticleData &a) const;
  inline bool operator==(int id_a) const;
  inline bool operator<(int id_a) const;

 private:
  /// Each particle has a unique identifier
  int id_;
  /// pdg id of the particle
  PdgCode pdgcode_;
  /// counter of the last collision/decay
  int id_process_;
  /// collision time
  double collision_time_;
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
inline void ParticleData::set_id(int i) {
  id_ = i;
}

/// look up the pdgcode of the particle
inline PdgCode ParticleData::pdgcode(void) const {
  return pdgcode_;
}

/// set id of the particle
inline void ParticleData::set_pdgcode(PdgCode i) {
  pdgcode_ = i;
}

/// look up the id of the collision process
inline int ParticleData::id_process(void) const {
  return id_process_;
}

/// set the id of the collision process
inline void ParticleData::set_id_process(int process_id) {
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
inline void ParticleData::set_collision_past(int id_counter) {
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

/// set particle four momentum by components
inline void ParticleData::set_momentum(const double &mass, const double &px,
                          const double &py, const double &pz) {
  momentum_.set_FourVector(sqrt(mass * mass + px * px + py * py + pz * pz),
                           px, py, pz);
}

/// particle position in Minkowski space
inline const FourVector &ParticleData::position(void) const {
  return position_;
}

/// set the particle position directly
inline void ParticleData::set_position(const FourVector &pos) {
  position_ = pos;
}

/// set the particle position by components
inline void ParticleData::set_position(const double &x0, const double &x1,
                          const double &x2, const double &x3) {
  position_.set_FourVector(x0, x1, x2, x3);
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
inline bool ParticleData::operator==(int id_a) const {
  return this->id_ == id_a;
}

/// sort particles along to given id
inline bool ParticleData::operator<(int id_a) const {
  return this->id_ < id_a;
}

}  // namespace Smash

#endif  // SRC_INCLUDE_PARTICLEDATA_H_
