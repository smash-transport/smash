/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PARTICLEDATA_H_
#define SRC_INCLUDE_PARTICLEDATA_H_

#include "FourVector.h"

class ParticleData {
  public:
  /* Use improbable values for default constructor */
  ParticleData() :id_(-1) {}
  int id(void);
  void inline set(const int &id, const double &momenta_l,
                  const double &momenta_t);
  void inline set_id(const int &id);
  double inline momenta_l(void);
  void inline set_momenta_l(const double &momenta_l);
  double inline momenta_t(void);
  void inline set_momenta_t(const double &momenta_t);
  void inline set_momenta(const double &momenta_l, const double &momenta_t);
  FourVector inline x(void);
  void inline set_position(const FourVector &position);

  private:
    /* Each particle has a unique identifier */
    int id_;
    /* longitudinal momenta */
    double momenta_l_;
    /* tranverse momenta */
    double momenta_t_;
    /* position in space */
    FourVector x_;
};

int inline ParticleData::id(void) {
  return id_;
}

void inline ParticleData::set_id(const int &i) {
  id_ = i;
}

double inline ParticleData::momenta_l(void) {
  return momenta_l_;
}

void inline ParticleData::set_momenta_l(const double &m_l) {
  momenta_l_ = m_l;
}

double inline ParticleData::momenta_t(void) {
  return momenta_t_;
}

void inline ParticleData::set_momenta_t(const double &m_t) {
  momenta_t_ = m_t;
}

void inline ParticleData::set_momenta(const double &m_l, const double &m_t) {
  momenta_l_ = m_l;
  momenta_t_ = m_t;
}

FourVector inline ParticleData::x(void) {
  return x_;
}

void inline ParticleData::set_position(const FourVector &pos) {
  x_ = pos;
}

#endif  // SRC_INCLUDE_PARTICLEDATA_H_
