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
  /* Explicit constructor */
  ParticleData(int i, double m_l, double m_t) : id_(i), momenta_l_(m_l),
    momenta_t_(m_t) {}
  int id(void);
  void inline set(int id, double momenta_l, double momenta_t);
  double inline momenta_l(void);
  void inline set_mometa_l(double momenta_l);
  double inline momenta_t(void);
  void inline set_mometa_t(double momenta_t);

  private:
    /* Each particle has a unique identifier */
    int id_;
    /* longitudinal momenta */
    double momenta_l_;
    /* tranverse momenta */
    double momenta_t_;
};

int inline ParticleData::id(void) {
  return id_;
}

void inline ParticleData::set(int i, double m_l, double m_t) {
  id_ = i;
  momenta_l_ = m_l;
  momenta_t_ = m_t;
}

double inline ParticleData::momenta_l(void) {
  return momenta_l_;
}

void inline ParticleData::set_mometa_l(double m_l) {
  momenta_l_ = m_l;
}

double inline ParticleData::momenta_t(void) {
  return momenta_t_;
}

void inline ParticleData::set_mometa_t(double m_t) {
  momenta_t_ = m_t;
}

#endif  // SRC_INCLUDE_PARTICLEDATA_H_
