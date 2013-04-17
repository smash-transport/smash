/*
 *    Copyright (c) 2012 maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PARTICLETYPE_H_
#define SRC_INCLUDE_PARTICLETYPE_H_

#include <string>

class ParticleType {
  public:
    /* Use improbable values for default constructor */
    ParticleType() : name_("unknown"), mass_(-1), lifetime_(-1),
      isospin_(100), pdgcode_(-1) {}
    /* Explicit constructor */
    ParticleType(std::string n, double m, int id) : name_(n), mass_(m),
      pdgcode_(id) {}
    /* access data */
    std::string inline name(void) const;
    double inline mass(void) const;
    int inline pdgcode(void) const;

  private:
    /* Data of the particle type */
    const std::string name_;
    const double mass_;
    double lifetime_;
    float isospin_;
    int pdgcode_;
};

double inline ParticleType::mass(void) const {
  return mass_;
}

std::string inline ParticleType::name(void) const {
  return name_;
}

int inline ParticleType::pdgcode(void) const {
  return pdgcode_;
}

#endif  // SRC_INCLUDE_PARTICLETYPE_H_
