/*
 *    Copyright (c) 2012 maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PARTICLETYPE_H_
#define SRC_INCLUDE_PARTICLETYPE_H_

#include <string>

class ParticleType {
  private:
    /* Data of the particle type */
    string name;
    double mass;
    float isospin;
};

#endif  // SRC_INCLUDE_PARTICLETYPE_H_
