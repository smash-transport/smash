/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PARTICLEDATA_H_
#define SRC_INCLUDE_PARTICLEDATA_H_

#include <math.h>

#include "../include/FourVector.h"

class ParticleData {
  public:
  /* Use improbable values for default constructor */
  ParticleData() :id_(-1), collision_id_(-1) {}
  void inline set(const int &id, const double &momenta_l,
                  const double &momenta_t);
  int id(void) const;
  void inline set_id(const int &id);
  int collision_id(void) const;
  void inline set_collision_id(const int &collision_id);
  double collision_time(void) const;
  void inline set_collision_time(const double &collision_time);
  void inline set_collision(const double &collision_time,
    const int &collision_id);
  FourVector inline momentum(void) const;
  void inline set_momentum(const FourVector &momentum_vector);
  void inline set_momentum(const double &mass, const double &px,
                           const double &py, const double &pz);
  FourVector inline x(void) const;
  void inline set_position(const FourVector &position);
  void inline set_position(const double &x0, const double &x3,
                           const double &x1, const double &x2);
  /* get velocities */
  double inline velocity_x(void) { return momentum().x2() / momentum().x0(); }
  double inline velocity_y(void) { return momentum().x3() / momentum().x0(); }
  double inline velocity_z(void) { return momentum().x1() / momentum().x0(); }
  /* overloaded operators */
  bool inline operator==(const ParticleData &a);

  private:
    /* Each particle has a unique identifier */
    int id_;
    /* Next particle we'd collide against */
    int collision_id_;
    /* collision time */
    double collision_time_;
    /* momenta of the particle */
    FourVector momentum_;
    /* position in space */
    FourVector x_;
};

int inline ParticleData::id(void) const {
  return id_;
}

void inline ParticleData::set_id(const int &i) {
  id_ = i;
}

int inline ParticleData::collision_id(void) const {
  return collision_id_;
}

void inline ParticleData::set_collision_id(const int &collision_i) {
  collision_id_ = collision_i;
}

double inline ParticleData::collision_time(void) const {
  return collision_time_;
}

void inline ParticleData::set_collision_time(const double &collision_t) {
  collision_time_ = collision_t;
}

void inline ParticleData::set_collision(const double &collision_t,
  const int &collision_i) {
  collision_time_ = collision_t;
  collision_id_ = collision_i;
}

FourVector inline ParticleData::momentum(void) const {
  return momentum_;
}

void inline ParticleData::set_momentum(const FourVector &momentum_vector) {
  momentum_ = momentum_vector;
}

void inline ParticleData::set_momentum(const double &mass, const double &px,
                          const double &py, const double &pz) {
  momentum_.set_FourVector(sqrt(mass * mass + px * px + py * py + pz * pz),
                           px, py, pz);
}

FourVector inline ParticleData::x(void) const {
  return x_;
}

void inline ParticleData::set_position(const FourVector &pos) {
  x_ = pos;
}

void inline ParticleData::set_position(const double &x0, const double &x3,
                          const double &x1, const double &x2) {
  x_.set_FourVector(x0, x3, x1, x2);
}

bool inline ParticleData::operator==(const ParticleData &a) {
  return this->id_ == a.id_;
}

#endif  // SRC_INCLUDE_PARTICLEDATA_H_
