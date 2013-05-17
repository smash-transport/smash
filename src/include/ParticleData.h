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
  ParticleData() :id_(-1), collision_id_(-1), collision_time_(0.0) {}
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
  FourVector inline position(void) const;
  void inline set_position(const FourVector &position);
  void inline set_position(const double &x0, const double &x1,
                           const double &x2, const double &x3);
  /* get velocities */
  double inline velocity_x(void) { return momentum().x1() / momentum().x0(); }
  double inline velocity_y(void) { return momentum().x2() / momentum().x0(); }
  double inline velocity_z(void) { return momentum().x3() / momentum().x0(); }
  /* overloaded operators */
  bool inline operator==(const ParticleData &a);

  private:
    /* Each particle has a unique identifier */
    int id_;
    /* Next particle we'd collide against */
    int collision_id_;
    /* collision time */
    double collision_time_;
    /* momenta of the particle: x0, x1, x2, x3 as E, px, py, pz */
    FourVector momentum_;
    /* position in space: x0, x1, x2, x3 as t, x, y, z */
    FourVector position_;
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

/* set particle four momentum directly */
void inline ParticleData::set_momentum(const FourVector &momentum_vector) {
  momentum_ = momentum_vector;
}

/* set particle four momentum by components */
void inline ParticleData::set_momentum(const double &mass, const double &px,
                          const double &py, const double &pz) {
  momentum_.set_FourVector(sqrt(mass * mass + px * px + py * py + pz * pz),
                           px, py, pz);
}

/* particle position in Minkowski space */
FourVector inline ParticleData::position(void) const {
  return position_;
}

/* set the particle position directly */
void inline ParticleData::set_position(const FourVector &pos) {
  position_ = pos;
}

/* set the particle position by components */
void inline ParticleData::set_position(const double &x0, const double &x3,
                          const double &x1, const double &x2) {
  position_.set_FourVector(x0, x3, x1, x2);
}

/* check if the particles are identical */
bool inline ParticleData::operator==(const ParticleData &a) {
  return this->id_ == a.id_;
}

#endif  // SRC_INCLUDE_PARTICLEDATA_H_
