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
#include "include/parameters.h"

struct ExperimentParameters;

class NucleusModus : public ModusDefault {
 public:
  /* default constructor with probable values */
  NucleusModus() = default;

  /* special class funtions */
  void assign_params(
      std::list<Parameters> *configuration);  // TODO(mkretz) -> ctor
  void print_startup();  // TODO(mkretz): needs to be discoverable from an
                         // outside "printer"

  void initial_conditions(Particles *particles,
                          const ExperimentParameters &parameters);

  float nuclear_radius(const int &A) const;
  float woods_saxon(const float &radius);
 
  // in ModusDefault:
  // * sanity_check
  // * check_collision_geometry
  // * propagate

 private:
  // Projectile and target particles PDG ID
  std::vector<int> projectile_;
  std::vector<int> target_;
  // Center-of-mass energy of the collision
  float sqrts_ = 1.f;
  // initial separation of nuclei from origin
  float initial_displacement_ = 1.1;
  // steepness of woods-saxon distribution
  float steepness = .545;
};

#endif  // SRC_INCLUDE_COLLIDERMODUS_H_
