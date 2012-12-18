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
  double momenta_t(void);
  double momenta_l(void);

  private:
    /* Each particle has a unique identifier */
    const int id_;
    /* longitudinal momenta */
    double momenta_l_;
    /* tranverse momenta */
    double momenta_t_;
};

#endif  // SRC_INCLUDE_PARTICLEDATA_H_
