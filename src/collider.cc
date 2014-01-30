/*
 *    Copyright (c) 2014
 *      Hannah Petersen <petersen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
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
    Modus::assign_params(configuration);
    bool match = false;
    std::list<Parameters>::iterator i = configuration->begin();
    while (i != configuration->end()) {
        char *key = i->key();
        char *value = i->value();
        printd("%s %s\n", key, value);
        /* string values */
        if (strcmp(key, "PROJECTILE") == 0) {
          strncpy(projectile_, value, sizeof(&projectile_));
            match = true;
        }
        if (strcmp(key, "TARGET") == 0) {
          strncpy(target_, value, sizeof(&target_));
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
    Modus::print_startup();
    printf("Projectile name: %s \n", projectile_);
    printf("Target name: %s \n", target_);
    printf("Center-of-mass energy %f GeV\n", sqrts_);
}

/* initial_conditions - sets particle data for @particles */
void ColliderModus::initial_conditions(Particles *particles) {
  /* velocity of particles */
  float beta;

  

}














