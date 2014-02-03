/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/Collider.h"
#include "include/CrossSections.h"
#include "include/Particles.h"
#include "include/constants.h"
#include "include/collisions.h"
#include "include/decays.h"
#include "include/distributions.h"
#include "include/input-decaymodes.h"
#include "include/input-particles.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/param-reader.h"
#include "include/Angles.h"

void ColliderModus::assign_params(std::list<Parameters>
                                          *configuration) {
    ModusDefault::assign_params(configuration);
    bool match = false;
    std::list<Parameters>::iterator i = configuration->begin();
    while (i != configuration->end()) {
        char *key = i->key();
        char *value = i->value();
        printd("%s %s\n", key, value);
        /* integer values */
        if (strcmp(key, "PROJECTILE") == 0) {
          projectile_ = abs(atoi(value));
            match = true;
        }
        if (strcmp(key, "TARGET") == 0) {
          target_ = abs(atoi(value));
            match = true;
        }
        /* float values */
        if (strcmp(key, "SQRTS") == 0) {
            sqrts_ = (fabs(atof(value)));
            match = true;
        }
        /* remove processed entry */
        if (match) {
            i = configuration->erase(i);
            match = false;
        } else {
            ++i;
        }
    }
}


/* print_startup - console output on startup of box specific parameters */
void ColliderModus::print_startup() {
    ModusDefault::print_startup();
    printf("Projectile PDG ID: %d \n", projectile_);
    printf("Target PDG ID: %d \n", target_);
    printf("Center-of-mass energy %10.3f GeV\n", sqrts_);
}

/* initial_conditions - sets particle data for @particles */
void ColliderModus::initial_conditions(Particles *particles) {
  particles->create(1, projectile_);
  particles->create(1, target_);

  for (auto i = particles->begin(); i !=particles->end(); i++) {
    float mass = particles->type(i->first).mass();
    printf("id %d pdgcode %d mass %f\n", i->first, i->second.pdgcode(), mass);
    /* velocity of particles */
    double cms_gamma = sqrts_ / mass;
    double cms_beta = sqrt(sqrts_ * sqrts_ - mass * mass / sqrts_);
    if (i->first == 0) {
     i->second.set_position(1.0, 0.0, 0.0, -1.0);
     i->second.set_momentum(mass, 0.0, 0.0, cms_gamma * cms_beta * mass);
    } else if (i->first == 1) {
     i->second.set_position(1.0, 0.0, 0.0, 1.0);
     i->second.set_momentum(mass, 0.0, 0.0, - cms_gamma *cms_beta * mass);
    }
  }
}
