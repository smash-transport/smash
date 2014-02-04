/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 * 
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include <cinttypes>

#include "include/ModusDefault.h"
#include "include/collisions.h"
#include "include/constants.h"
#include "include/Experiment.h"
#include "include/outputroutines.h"

/*general propagation routine */
void ModusDefault::propagate(Particles *particles, const ExperimentParameters &parameters) {
    FourVector distance, position;
    for (auto i = particles->begin(); i != particles->end(); ++i) {
        /* propagation for this time step */
        distance.set_FourVector(parameters.eps,
                                i->second.velocity_x() * parameters.eps,
                                i->second.velocity_y() * parameters.eps,
                                i->second.velocity_z() * parameters.eps);
        printd("Particle %d motion: %g %g %g %g\n", i->first,
               distance.x0(), distance.x1(), distance.x2(), distance.x3());
        position = i->second.position();
        position += distance;
        i->second.set_position(position);
    }
}

// check particle pairs for collision
void ModusDefault::check_collision_geometry(
    Particles *particles, CrossSections *cross_sections,
    std::list<int> *collision_list, size_t *rejection_conflict,
    const ExperimentParameters &parameters) {
    FourVector distance;
    double radial_interaction = sqrt(parameters.cross_section * fm2_mb
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
                                        collision_list, parameters.eps,
                                        i->first, j->first, rejection_conflict);
        }
    }
}
