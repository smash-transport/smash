/*
 *
 *    Copyright (c) 2012-2013
 *      Hannah Petersen <petersen@fias.uni-frankfurt.de>
 *      
 *    GNU General Public License (GPLv3)
 *
 */

#include <map>
#include <list>

#include "include/Box.h"
#include "include/Experiment.h"
#include "include/CrossSections.h"
#include "include/Modus.h"
#include "include/Parameters.h"
#include "include/collisions.h"
#include "include/decays.h"
#include "include/input-particles.h"
#include "include/input-decaymodes.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/param-reader.h"

/* #include "include/Sphere.h" */

/* Experiment carries everything that is needed for the evolution */
std::unique_ptr<Experiment> Experiment::create(char *modus_chooser) {
  typedef std::unique_ptr<Experiment> ExperimentPointer;
  if (strcmp(modus_chooser, "Box") == 0) {
    return ExperimentPointer{new ExperimentImplementation<BoxModus>};
  }
//  else if (modus == 2) {
//    return ExperimentPointer{new ExperimentImplementation<SphereModus>};
//  }
    else {
    throw std::string("Invalid Modus requested from Experiment::create.");
  }
}

/*This method reads the parameters in */
template <typename Modus>
void ExperimentImplementation<Modus>::configure(std::list<Parameters>
                                                configuration) {
    bc.assign_params(&configuration);
    warn_wrong_params(&configuration);
    bc.print_startup();
/* reducing cross section according to number of test particle */
    if (bc.testparticles > 1) {
      printf("IC test particle: %i\n", bc.testparticles);
      bc.cross_section = bc.cross_section / bc.testparticles;
      printf("Elastic cross section: %g mb\n", bc.cross_section);
    }
}

/* This method reads the particle type and cross section information
 * and does the initialization of the system (fill the particles map)*/
template <typename Modus>
void ExperimentImplementation<Modus>::initialize(char *path) {
    srand48(bc.seed);
    input_particles(particles, path);
    input_decaymodes(particles, path);
    cross_sections->add_elastic_parameter(bc.cross_section);
    bc.initial_conditions(particles);
    bc.energy_initial = bc.energy_total(particles);
    write_measurements_header(*particles);
    print_header();
    write_particles(*particles);
}

/*This is the loop over timesteps, carrying out collisions and decays 
 * and propagating particles */
template <typename Modus>
void ExperimentImplementation<Modus>::run_time_evolution() {
    bc.sanity_check(particles);
    std::list<int> collision_list, decay_list;
    size_t interactions_total = 0, previous_interactions_total = 0,
    interactions_this_interval = 0;
    size_t rejection_conflict = 0;
    int resonances = 0, decays = 0;
    print_measurements(*particles, interactions_total,
                       interactions_this_interval, bc.energy_initial,
                       bc.time_start);
    for (int step = 0; step < bc.steps; step++) {
        /* Check resonances for decays */
        check_decays(particles, &decay_list, bc);
        /* Do the decays */
        if (!decay_list.empty()) {
            decays += decay_list.size();
            interactions_total = decay_particles(particles,
                                       &decay_list, interactions_total);
        }
        /* fill collision table by cells */
        bc.check_collision_geometry(particles, cross_sections,
                                 &collision_list, &rejection_conflict);
        /* particle interactions */
        if (!collision_list.empty()) {
            printd_list(collision_list);
            interactions_total = collide_particles(particles, &collision_list,
                                             interactions_total, &resonances);
        }
        bc.propagate(particles);
        /* physics output during the run */
        if (step > 0 && (step + 1) % bc.output_interval == 0) {
            interactions_this_interval = interactions_total
            - previous_interactions_total;
            previous_interactions_total = interactions_total;
            print_measurements(*particles, interactions_total,
                              interactions_this_interval, bc.energy_initial,
                              bc.time_start);
            printd("Resonances: %i Decays: %i\n", resonances, decays);
            printd("Ignored collisions %zu\n", rejection_conflict);
            /* save evolution data */
            write_measurements(*particles, interactions_total,
                               interactions_this_interval, resonances, decays,
                               rejection_conflict);
            write_vtk(*particles);
        }
    }
        /* Guard against evolution */
        if (likely(bc.steps > 0)) {
        /* if there are no particles no interactions happened */
            if (likely(!particles->empty())) {
             print_tail(bc.time_start, interactions_total * 2
                        / particles->time() / particles->size());
            } else {
             print_tail(bc.time_start, 0);
             printf("Total ignored collisions: %zu\n", rejection_conflict);
            }
        }
}

/* Tear down of everything */
template <typename Modus>
void ExperimentImplementation<Modus>::end() {
    delete particles;
    delete cross_sections;
}

