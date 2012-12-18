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
    ParticleType() : name_("unknown"), mass_(-1) {}
    /* Explicit constructor */
    ParticleType(std::string n, double m) : name_(n), mass_(m) {}
    /* access data */
    std::string inline name(void);
    double inline mass(void);

  private:
    /* Data of the particle type */
    const std::string name_;
    const double mass_;
    double lifetime_;
    float isospin_;
};

double inline ParticleType::mass(void) {
  return mass_;
}

std::string inline ParticleType::name(void) {
  return name_;
}

#endif  // SRC_INCLUDE_PARTICLETYPE_H_
