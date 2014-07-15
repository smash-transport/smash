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

#include "forwarddeclarations.h"
#include "modusdefault.h"

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
    
  float initial_conditions(Particles *particles,
                          const ExperimentParameters &parameters);
  
 private:
  /* Sphere radius length */
  float radius_;
  /* Total number of particles in Sphere */
  int number_of_particles_;
  /* Temperature for momentum distribution */
  float sphere_temperature_;
  /// Starting time for the Sphere
    const float start_time_ = 0.0f;    
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SPHEREMODUS_H_
