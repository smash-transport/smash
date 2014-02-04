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

#include "../include/CrossSections.h"
#include "../include/ModusDefault.h"
#include "../include/Particles.h"
#include "../include/Parameters.h"

class ExperimentParameters;

class ColliderModus : public ModusDefault {
 public:
  /* default constructor with probable values */
  ColliderModus() = default;

    /* special class funtions */
    // XXX: -> ctor
    void assign_params(std::list<Parameters> *configuration);
    // XXX: needs to be discoverable from an outside "printer"
    void print_startup();

    void initial_conditions(Particles *particles,
                            const ExperimentParameters &parameters);

    // in ModusDefault:
    // * sanity_check
    // * check_collision_geometry
    // * propagate

 private:
    /* Projectile particle PDG ID*/
    int projectile_ = 2212;
    /* Target particle PDG ID*/
    int target_ = 2212;
    /* Center-of-mass energy of the collision */
    float sqrts_ = 1.f;
};

#endif  // SRC_INCLUDE_COLLIDERMODUS_H_
