/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_COLLIDERMODUS_H_
#define SRC_INCLUDE_COLLIDERMODUS_H_

#include <stdint.h>
#include <cmath>
#include <list>
#include <string>

#include "include/crosssections.h"
#include "include/modusdefault.h"
#include "include/particles.h"
#include "include/random.h"

namespace Smash {

class Configuration;
struct ExperimentParameters;

class ColliderModus : public ModusDefault {
 public:
  /* default constructor with probable values */
  explicit ColliderModus(Configuration modus_config,
           const ExperimentParameters &parameters);

  void print_startup();  // TODO(mkretz): needs to be discoverable from an
                         // outside "printer"

  float initial_conditions(Particles *particles,
                          const ExperimentParameters &parameters);

  // in ModusDefault:
  // * sanity_check
  // * check_collision_geometry
  // * propagate

 private:
  /* Projectile particle PDG ID*/
  //TODO(mkretz): Matthias wants to fix this back to const.
  PdgCode projectile_;
  /* Target particle PDG ID*/
  PdgCode target_;
  /* Center-of-mass energy of the collision */
  const float sqrts_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_COLLIDERMODUS_H_
