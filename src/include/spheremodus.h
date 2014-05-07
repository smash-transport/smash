/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SPHEREMODUS_H_
#define SRC_INCLUDE_SPHEREMODUS_H_

#include <stdint.h>
#include <cmath>
#include <list>

#include "include/crosssections.h"
#include "include/modusdefault.h"
#include "include/particles.h"

namespace Smash {

class SphereModus;
class Configuration;
struct ExperimentParameters;
    
class SphereModus : public ModusDefault {
 public:
  /* default constructor with probable values */
  explicit SphereModus(Configuration modus_config,
           const ExperimentParameters &parameters);
  
  void print_startup();
    
  void initial_conditions(Particles *particles,
                          const ExperimentParameters &parameters);
  
 private:
  /* Sphere radius length */
  float radius_;
  /* Total number of particles in Sphere */
  int number_of_particles_; 
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SPHEREMODUS_H_
