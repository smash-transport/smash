/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include <cmath>
#include <cstdlib>
#include <list>

#include "include/collidermodus.h"
#include "include/configuration.h"
#include "include/experimentparameters.h"
#include "include/outputroutines.h"
#include "include/parameters.h"

ColliderModus::ColliderModus(Configuration &config)
    : projectile_(config.take({"Collider", "PROJECTILE"})),
      target_    (config.take({"Collider", "TARGET"})),
      sqrts_     (config.take({"Collider", "SQRTS"})) {
}

/* print_startup - console output on startup of box specific parameters */
void ColliderModus::print_startup() {
  printf("Projectile PDG ID: %d \n", projectile_);
  printf("Target PDG ID: %d \n", target_);
  printf("Center-of-mass energy %10.3f GeV\n", sqrts_);
}

/* initial_conditions - sets particle data for @particles */
void ColliderModus::initial_conditions(Particles *particles,
                                       const ExperimentParameters &) {
  particles->create(1, projectile_);
  particles->create(1, target_);

  for (auto i = particles->begin(); i != particles->end(); i++) {
    float mass = particles->type(i->first).mass();
    printf("id %d pdgcode %d mass %f\n", i->first, i->second.pdgcode(), mass);
    /* velocity of particles */
    double cms_gamma = sqrts_ / mass;
    double cms_beta = sqrt(sqrts_ * sqrts_ - mass * mass / sqrts_);
    // Sample impact parameter
    double impact_parameter = drand48() * 5.0;
    if (i->first == 0) {
      i->second.set_position(1.0, impact_parameter, 0.0, -1.0);
      i->second.set_momentum(mass, 0.0, 0.0, cms_gamma * cms_beta * mass);
    } else if (i->first == 1) {
      i->second.set_position(1.0, 0.0, 0.0, 1.0);
      i->second.set_momentum(mass, 0.0, 0.0, -cms_gamma * cms_beta * mass);
    }
  }
}
