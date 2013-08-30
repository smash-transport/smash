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
  ParticleData(int i) :id_(i), pdgcode_(-1), id_partner_(-1), id_process_(-1),
    collision_time_(0.0), process_type_(-1) {}
  void inline set(int id, const double &momenta_l, const double &momenta_t);
  int id(void) const;
  void inline set_id(int id);
  int pdgcode(void) const;
  void inline set_pdgcode(int pdgcode);
  int id_partner(void) const;
  void inline set_id_partner(int id_b);
  int id_process(void) const;
  void inline set_id_process(int id);
  double collision_time(void) const;
  int process_type(void) const;
  void inline set_collision_time(const double &collision_time);
  void inline set_collision(int process_type, const double &collision_time,
                            int collision_id);
  void inline set_collision_past(int process_id);
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
  bool inline operator==(const ParticleData &a) const;
  bool inline operator<(const ParticleData &a) const;
  bool inline operator==(int id_a) const;
  bool inline operator<(int id_a) const;

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
int inline ParticleData::id(void) const {
  return id_;
}

/* set id of the particle */
void inline ParticleData::set_id(int i) {
  id_ = i;
}

/* look up the pdgcode of the particle */
int inline ParticleData::pdgcode(void) const {
  return pdgcode_;
}

/* set id of the particle */
void inline ParticleData::set_pdgcode(int i) {
  pdgcode_ = i;
}

/* look up the id of the collision partner */
int inline ParticleData::id_partner(void) const {
  return id_partner_;
}

/* set the id of the collision partner */
void inline ParticleData::set_id_partner(int id_b) {
  id_partner_ = id_b;
}

/* look up the id of the collision process */
int inline ParticleData::id_process(void) const {
  return id_process_;
}

/* set the id of the collision process */
void inline ParticleData::set_id_process(int process_id) {
  id_process_ = process_id;
}

double inline ParticleData::collision_time(void) const {
  return collision_time_;
}

int inline ParticleData::process_type(void) const {
  return process_type_;
}

void inline ParticleData::set_collision_time(const double &collision_t) {
  collision_time_ = collision_t;
}

/* set possible collision data */
void inline ParticleData::set_collision(int proc_type,
  const double &collision_t, int id_b) {
  process_type_ = proc_type;
  collision_time_ = collision_t;
  id_partner_ = id_b;
}

/* set happened collision data */
void inline ParticleData::set_collision_past(int id_counter) {
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
bool inline ParticleData::operator==(const ParticleData &a) const {
  return this->id_ == a.id_;
}

/* sort particles along their id */
bool inline ParticleData::operator<(const ParticleData &a) const {
  return this->id_ < a.id_;
}

/* check if the particles are identical to a given id */
bool inline ParticleData::operator==(int id_a) const {
  return this->id_ == id_a;
}

/* sort particles along to given id */
bool inline ParticleData::operator<(int id_a) const {
  return this->id_ < id_a;
}

#endif  // SRC_INCLUDE_PARTICLEDATA_H_
