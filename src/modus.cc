/*
 *
 *    Copyright (c) 2013-2014
 *      Hannah Petersen <petersen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include <cinttypes>

#include "include/Modus.h"
#include "include/collisions.h"
#include "include/constants.h"
#include "include/outputroutines.h"

void Modus::assign_params(std::list<Parameters> *configuration) {
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
void Modus::print_startup() {
    printf("Elastic cross section: %g mb\n", cross_section);
    printf("Using temporal stepsize: %g fm/c\n", eps);
    printf("Maximum number of steps: %i \n", steps);
    printf("Random number seed: %" PRId64 "\n", seed);
}

/* calculates the total energy in the system from zero component of
 * all momenta of particles
 * XXX should be expanded to all quantum numbers of interest */
float Modus::energy_total(Particles *particles) {
    float energy_sum = 0.0;
    for (auto i = particles->begin(); i != particles->end(); ++i) {
         energy_sum += i->second.momentum().x0();
     }
    return energy_sum;
}

/*general propagation routine */
void Modus::propagate(Particles *particles) {
    FourVector distance, position;
    for (auto i = particles->begin(); i != particles->end(); ++i) {
        /* propagation for this time step */
        distance.set_FourVector(eps,
                                i->second.velocity_x() * eps,
                                i->second.velocity_y() * eps,
                                i->second.velocity_z() * eps);
        printd("Particle %d motion: %g %g %g %g\n", i->first,
               distance.x0(), distance.x1(), distance.x2(), distance.x3());
    }
}

/*empty methods are needed in Boxmodus */
int Modus::sanity_check(Particles *particles __attribute__((unused))) {
    return 0;
}

FourVector Modus::boundary_condition(FourVector position,
                  bool *boundary_hit __attribute__((unused))) {
    return position;
}

// check particle pairs for collision
void Modus::check_collision_geometry(Particles *particles,
       CrossSections *cross_sections, std::list<int> *collision_list,
       size_t *rejection_conflict) {
    FourVector distance;
    double radial_interaction = sqrt(cross_section * fm2_mb
                                     * M_1_PI) * 2;
    for (auto i = particles->begin(); i != particles->end(); ++i) {
        for (auto j = particles->begin(); j != particles->end(); ++j) {
            /* exclude check on same particle and double counting */
            if (i->first >= j->first)
                continue;
            distance = i->second.position() - j->second.position();
            /* skip particles that are double interaction radius length away */
            if (distance > radial_interaction)
                continue;
            collision_criteria_geometry(particles, cross_sections,
                                        collision_list, this->eps,
                                        i->first, j->first, rejection_conflict);
        }
    }
}
