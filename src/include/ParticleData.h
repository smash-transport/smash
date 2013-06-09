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
  ParticleData() :id_(-1), id_partner_(-1), id_process_(-1),
    collision_time_(0.0), process_type_(-1), lifetime_(-1) {}
  void inline set(const int &id, const double &momenta_l,
                  const double &momenta_t);
  int id(void) const;
  void inline set_id(const int &id);
  int id_partner(void) const;
  void inline set_id_partner(const int &id_b);
  int id_process(void) const;
  void inline set_id_process(const int &id);
  double collision_time(void) const;
  int process_type(void) const;
  double lifetime(void) const;
  void inline set_lifetime(const double &time);
  void inline set_collision_time(const double &collision_time);
  void inline set_collision(const int &process_type,
                            const double &collision_time,
                            const int &collision_id);
  void inline set_collision_past(const int &process_id);
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
    int id_partner_;
    /* counter of the last collision/decay */
    int id_process_;
    /* collision time */
    double collision_time_;
    /* Type of interaction. 0: 2->2, 1: 2->1, 2: 1->2 */
    int process_type_;
    /* momenta of the particle: x0, x1, x2, x3 as E, px, py, pz */
    FourVector momentum_;
    /* position in space: x0, x1, x2, x3 as t, x, y, z */
    FourVector position_;
    /* particle lifetime */
    double lifetime_;
};

int inline ParticleData::id(void) const {
  return id_;
}

void inline ParticleData::set_id(const int &i) {
  id_ = i;
}

/* look up the id of the collision partner */
int inline ParticleData::id_partner(void) const {
  return id_partner_;
}

/* set the id of the collision partner */
void inline ParticleData::set_id_partner(const int &id_b) {
  id_partner_ = id_b;
}

/* look up the id of the collision process */
int inline ParticleData::id_process(void) const {
  return id_process_;
}

/* set the id of the collision process */
void inline ParticleData::set_id_process(const int &process_id) {
  id_process_ = process_id;
}

double inline ParticleData::collision_time(void) const {
  return collision_time_;
}

int inline ParticleData::process_type(void) const {
  return process_type_;
}

double inline ParticleData::lifetime(void) const {
  return lifetime_;
}

void inline ParticleData::set_lifetime(const double &time) {
  lifetime_ = time;
}

void inline ParticleData::set_collision_time(const double &collision_t) {
  collision_time_ = collision_t;
}

/* set possible collision data */
void inline ParticleData::set_collision(const int &proc_type,
  const double &collision_t, const int &id_b) {
  process_type_ = proc_type;
  collision_time_ = collision_t;
  id_partner_ = id_b;
}

/* set happened collision data */
void inline ParticleData::set_collision_past(const int &id_counter) {
  collision_time_ = 0.0;
  id_process_ = id_counter;
  id_partner_ = -1;
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
