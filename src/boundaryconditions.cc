/*
 *
 *    Copyright (c) 2013-2014
 *      Hannah Petersen <petersen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */

#include "include/BoundaryConditions.h"
#include "include/outputroutines.h"

void BoundaryConditions::assign_params(std::list<Parameters> *configuration) {
    bool match = false;
    std::list<Parameters>::iterator i = configuration->begin();
    while (i != configuration->end()) {
        char *key = i->key();
        char *value = i->value();
        printd("Looking for match %s %s\n", key, value);
        /* integer values */
        if (strcmp(key, "STEPS") == 0) {
            steps = (abs(atoi(value)));
            match = true;
        }
        if (strcmp(key, "RANDOMSEED") == 0) {
            /* negative seed means random startup value */
            if (atol(value) > 0)
                seed = (atol(value));
            else
                seed = (time(NULL));
            match = true;
        }
        if (strcmp(key, "UPDATE") == 0) {
            output_interval = (abs(atoi(value)));
            match = true;
        }
        if (strcmp(key, "TESTPARTICLES") == 0) {
            testparticles = (abs(atoi(value)));
            match = true;
        }
        /* double or float values */
        if (strcmp(key, "EPS") == 0) {
            eps = (fabs(atof(value)));
            match = true;
        }
        if (strcmp(key, "SIGMA") == 0) {
            cross_section = (fabs(atof(value)));
            match = true;
        }
        /* remove processed entry */
        if (match) {
            printd("Erasing %s %s\n", key, value);
            i = configuration->erase(i);
            match = false;
        } else {
            ++i;
        }
    }
}

/* print_startup - console output on startup of general parameters */
void BoundaryConditions::print_startup() {
    printf("Elastic cross section: %g mb\n", cross_section);
    printf("Using temporal stepsize: %g fm/c\n", eps);
    printf("Maximum number of steps: %i \n", steps);
    printf("Random number seed: %lli \n", seed);
}

float BoundaryConditions::energy_total(Particles *particles) {
    float energy_total = 0.0;
    for (std::map<int, ParticleData>::iterator i = particles->begin();
         i != particles->end(); ++i) {
         energy_total += i->second.momentum().x0();
     }
    return energy_total;
}

void BoundaryConditions::propagate(Particles *particles) {
    FourVector distance, position;
    for (std::map<int, ParticleData>::iterator i = particles->begin();
         i != particles->end(); ++i) {
        /* propagation for this time step */
        distance.set_FourVector(eps,
                                i->second.velocity_x() * eps,
                                i->second.velocity_y() * eps,
                                i->second.velocity_z() * eps);
        printd("Particle %d motion: %g %g %g %g\n", i->first,
               distance.x0(), distance.x1(), distance.x2(), distance.x3());
    }
}

int BoundaryConditions::prepare_evolution(Particles *particles) {
    return 0;
}

FourVector BoundaryConditions::boundary_condition(FourVector position,
                                                  bool *boundary_hit) {
}

// XXX needs to be implemented in general form
void BoundaryConditions::check_collision_geometry(Particles *particles,
       CrossSections *cross_sections, std::list<int> *collision_list,
       size_t *rejection_conflict) {
}
