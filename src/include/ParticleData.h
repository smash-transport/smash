/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PARTICLEDATA_H_
#define SRC_INCLUDE_PARTICLEDATA_H_

#include "FourVector.h"

#include "math.h"

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

#endif  // SRC_INCLUDE_PARTICLEDATA_H_
