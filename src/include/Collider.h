/*
 *    Copyright (c) 2014
 *      Hannah Petersen <petersen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_COLLIDER_H_
#define SRC_INCLUDE_COLLIDER_H_

#include <stdint.h>
#include <cmath>
#include <list>

#include "../include/CrossSections.h"
#include "../include/Modus.h"
#include "../include/Particles.h"
#include "../include/Parameters.h"


class ColliderModus : public Modus {
 public:
  /* default constructor with probable values */
 ColliderModus(): projectile_("unknown"), target_("unknown"), sqrts_(1.0f) {}
    /* special class funtions */
    void assign_params(std::list<Parameters> *configuration);
    void print_startup();
    void initial_conditions(Particles *particles);
    int sanity_check(Particles *particles);
    void propagate(Particles *particles);
    FourVector boundary_condition(FourVector position,
                                  bool *boundary_hit);
    void check_collision_geometry(Particles *particles,
                  CrossSections *cross_sections, std::list<int> *collision_list,
                  size_t *rejection_conflict);

 private:
    /* Projectile particle */
    std:string projectile_;
    /* Target particle */
    std:string target_;
    /* Center-of-mass energy of the collision */
    float sqrts_;
};



#endif  // SRC_INCLUDE_COLLIDER_H_
