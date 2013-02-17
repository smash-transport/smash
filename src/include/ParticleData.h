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
  ParticleData() :id_(-1) {}
  int id(void);
  void inline set(const int &id, const double &momenta_l,
                  const double &momenta_t);
  void inline set_id(const int &id);
  FourVector inline momentum(void);
  void inline set_momentum(const FourVector &momentum_vector);
  void inline set_momentum(const double &mass, const double &px,
                           const double &py, const double &pz);
  FourVector inline x(void);
  void inline set_position(const FourVector &position);
  void inline set_position(const double &x0, const double &x3,
                           const double &x1, const double &x2);
  void inline add_position(const FourVector &position);
  /* get velocities */
  double inline velocity_x(void) { return momentum().x2() / momentum().x0(); }
  double inline velocity_y(void) { return momentum().x3() / momentum().x0(); }
  double inline velocity_z(void) { return momentum().x1() / momentum().x0(); }

  private:
    /* Each particle has a unique identifier */
    int id_;
    /* momenta of the particle */
    FourVector momentum_;
    /* position in space */
    FourVector x_;
};

int inline ParticleData::id(void) {
  return id_;
}

void inline ParticleData::set_id(const int &i) {
  id_ = i;
}

FourVector inline ParticleData::momentum(void) {
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

FourVector inline ParticleData::x(void) {
  return x_;
}

void inline ParticleData::set_position(const FourVector &pos) {
  x_ = pos;
}

void inline ParticleData::set_position(const double &x0, const double &x3,
                          const double &x1, const double &x2) {
  x_.set_FourVector(x0, x3, x1, x2);
}

void inline ParticleData::add_position(const FourVector &pos) {
  x_ += pos;
}

#endif  // SRC_INCLUDE_PARTICLEDATA_H_
