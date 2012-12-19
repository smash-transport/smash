/*
 *    Copyright (c) 2012 maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PARTICLEDATA_H_
#define SRC_INCLUDE_PARTICLEDATA_H_

class ParticleData {
  public:
  /* Use improbable values for default constructor */
  ParticleData() :id_(-1) {}
  int id(void);
  void inline set(int id, double momenta_l, double momenta_t);
  void inline set_id(int id);
  double inline momenta_l(void);
  void inline set_momenta_l(double momenta_l);
  double inline momenta_t(void);
  void inline set_momenta_t(double momenta_t);
  void inline set_momenta(double momenta_l, double momenta_t);

  private:
    /* Each particle has a unique identifier */
    int id_;
    /* longitudinal momenta */
    double momenta_l_;
    /* tranverse momenta */
    double momenta_t_;
    /* position in space */
    float x_;
    float y_;
    float z_;
};

int inline ParticleData::id(void) {
  return id_;
}

void inline ParticleData::set_id(int i) {
  id_ = i;
}

double inline ParticleData::momenta_l(void) {
  return momenta_l_;
}

void inline ParticleData::set_momenta_l(double m_l) {
  momenta_l_ = m_l;
}

double inline ParticleData::momenta_t(void) {
  return momenta_t_;
}

void inline ParticleData::set_momenta_t(double m_t) {
  momenta_t_ = m_t;
}

void inline ParticleData::set_momenta(double m_l, double m_t) {
  momenta_l_ = m_l;
  momenta_t_ = m_t;
}

#endif  // SRC_INCLUDE_PARTICLEDATA_H_
