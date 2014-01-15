/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PARTICLEDATA_H_
#define SRC_INCLUDE_PARTICLEDATA_H_

#include <math.h>

#include "../include/FourVector.h"

class ParticleData {
 public:
  /* Use improbable values for default constructor */
  ParticleData() :id_(-1), pdgcode_(-1), id_partner_(-1), id_process_(-1),
    collision_time_(0.0), process_type_(-1) {}
  explicit ParticleData(int i) :id_(i), pdgcode_(-1), id_partner_(-1),
    id_process_(-1), collision_time_(0.0), process_type_(-1) {}
  inline void set(int id, const double &momenta_l, const double &momenta_t);
  inline int id(void) const;
  inline void set_id(int id);
  inline int pdgcode(void) const;
  inline void set_pdgcode(int pdgcode);
  inline int id_partner(void) const;
  inline void set_id_partner(int id_b);
  inline int id_process(void) const;
  inline void set_id_process(int id);
  inline double collision_time(void) const;
  inline int process_type(void) const;
  inline void set_collision_time(const double &collision_time);
  inline void set_collision(int process_type, const double &collision_time,
                            int collision_id);
  inline void set_collision_past(int process_id);
  inline FourVector momentum(void) const;
  inline void set_momentum(const FourVector &momentum_vector);
  inline void set_momentum(const double &mass, const double &px,
                           const double &py, const double &pz);
  inline FourVector position(void) const;
  inline void set_position(const FourVector &position);
  inline void set_position(const double &x0, const double &x1,
                           const double &x2, const double &x3);
  /* get velocities */
  inline double velocity_x(void) { return momentum().x1() / momentum().x0(); }
  inline double velocity_y(void) { return momentum().x2() / momentum().x0(); }
  inline double velocity_z(void) { return momentum().x3() / momentum().x0(); }
  /* overloaded operators */
  inline bool operator==(const ParticleData &a) const;
  inline bool operator<(const ParticleData &a) const;
  inline bool operator==(int id_a) const;
  inline bool operator<(int id_a) const;

 private:
    /* Each particle has a unique identifier */
    int id_;
    /* pdg id of the particle */
    int pdgcode_;
    /* Next particle we'd collide against */
    int id_partner_;
    /* counter of the last collision/decay */
    int id_process_;
    /* collision time */
    double collision_time_;
    /* Type of interaction. 0: 2->2, 2: 1->2 >99 ( = PDG code): 2->1 */
    int process_type_;
    /* momenta of the particle: x0, x1, x2, x3 as E, px, py, pz */
    FourVector momentum_;
    /* position in space: x0, x1, x2, x3 as t, x, y, z */
    FourVector position_;
};

/* look up the id of the particle */
inline int ParticleData::id(void) const {
  return id_;
}

/* set id of the particle */
inline void ParticleData::set_id(int i) {
  id_ = i;
}

/* look up the pdgcode of the particle */
inline int ParticleData::pdgcode(void) const {
  return pdgcode_;
}

/* set id of the particle */
inline void ParticleData::set_pdgcode(int i) {
  pdgcode_ = i;
}

/* look up the id of the collision partner */
inline int ParticleData::id_partner(void) const {
  return id_partner_;
}

/* set the id of the collision partner */
inline void ParticleData::set_id_partner(int id_b) {
  id_partner_ = id_b;
}

/* look up the id of the collision process */
inline int ParticleData::id_process(void) const {
  return id_process_;
}

/* set the id of the collision process */
inline void ParticleData::set_id_process(int process_id) {
  id_process_ = process_id;
}

inline double ParticleData::collision_time(void) const {
  return collision_time_;
}

inline int ParticleData::process_type(void) const {
  return process_type_;
}

inline void ParticleData::set_collision_time(const double &collision_t) {
  collision_time_ = collision_t;
}

/* set possible collision data */
inline void ParticleData::set_collision(int proc_type,
  const double &collision_t, int id_b) {
  process_type_ = proc_type;
  collision_time_ = collision_t;
  id_partner_ = id_b;
}

/* set happened collision data */
inline void ParticleData::set_collision_past(int id_counter) {
  collision_time_ = 0.0;
  id_process_ = id_counter;
  id_partner_ = -1;
}

inline FourVector ParticleData::momentum(void) const {
  return momentum_;
}

/* set particle four momentum directly */
inline void ParticleData::set_momentum(const FourVector &momentum_vector) {
  momentum_ = momentum_vector;
}

/* set particle four momentum by components */
inline void ParticleData::set_momentum(const double &mass, const double &px,
                          const double &py, const double &pz) {
  momentum_.set_FourVector(sqrt(mass * mass + px * px + py * py + pz * pz),
                           px, py, pz);
}

/* particle position in Minkowski space */
inline FourVector ParticleData::position(void) const {
  return position_;
}

/* set the particle position directly */
inline void ParticleData::set_position(const FourVector &pos) {
  position_ = pos;
}

/* set the particle position by components */
inline void ParticleData::set_position(const double &x0, const double &x1,
                          const double &x2, const double &x3) {
  position_.set_FourVector(x0, x1, x2, x3);
}

/* check if the particles are identical */
inline bool ParticleData::operator==(const ParticleData &a) const {
  return this->id_ == a.id_;
}

/* sort particles along their id */
inline bool ParticleData::operator<(const ParticleData &a) const {
  return this->id_ < a.id_;
}

/* check if the particles are identical to a given id */
inline bool ParticleData::operator==(int id_a) const {
  return this->id_ == id_a;
}

/* sort particles along to given id */
inline bool ParticleData::operator<(int id_a) const {
  return this->id_ < id_a;
}

#endif  // SRC_INCLUDE_PARTICLEDATA_H_

