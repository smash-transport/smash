/*
 *
 *    Copyright (c) 2015-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *    If Pythia cite 
 *    T. Sjöstrand, S. Mrenna and P. Skands, JHEP05 (2006) 026,
 *                          Comput. Phys. Comm. 178 (2008) 852.
 *
 */

#ifndef SRC_INCLUDE_PYTHIA_H_
#define SRC_INCLUDE_PYTHIA_H_

#include "forwarddeclarations.h"

namespace Smash {
  bool sortfunc (ParticleList p1, ParticleList p2);
  ParticleList string_excitation(const ParticleList &incoming_particles_,
								 const float formation_time_);
}

#endif  // SRC_INCLUDE_PYTHIA_H_
