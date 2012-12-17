/*
 *    Copyright (c) 2012 maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PARTICLETYPE_H_
#define SRC_INCLUDE_PARTICLETYPE_H_

#include <string>

class ParticleType : public ParticleData {
  public:
    /* Use improbable values for default constructor */
    ParticleType() : name("unknown"), mass(-1) {}
    /* Explicit constructor */
    ParticleType(std::string n, double m) : name(n), mass(m) {}

  private:
    /* Data of the particle type */
    const std::string name;
    const double mass;
    double lifetime;
    float isospin;
};

#endif  // SRC_INCLUDE_PARTICLETYPE_H_
