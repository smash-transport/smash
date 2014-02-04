/*
 *
 *    Copyright (c) 2012-2013
 *      SMASH Team
 * 
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <map>
#include <list>
#include <cinttypes>

#include "include/BoxModus.h"
#include "include/Experiment.h"
#include "include/ColliderModus.h"
#include "include/CrossSections.h"
#include "include/ModusDefault.h"
#include "include/Parameters.h"
#include "include/collisions.h"
#include "include/decays.h"
#include "include/input-particles.h"
#include "include/input-decaymodes.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/param-reader.h"
#include "include/time.h"

/* #include "include/SphereModus.h" */

/* ExperimentBase carries everything that is needed for the evolution */
std::unique_ptr<ExperimentBase> ExperimentBase::create(std::string modus_chooser, int nevents) {
  typedef std::unique_ptr<ExperimentBase> ExperimentPointer;
  if (modus_chooser.compare("Box") == 0) {
    return ExperimentPointer(new Experiment<BoxModus>(nevents));
  } else if (modus_chooser.compare("Collider") == 0) {
    return ExperimentPointer(new Experiment<ColliderModus>(nevents));
  } else {
    throw std::string("Invalid ModusDefault requested from ExperimentBase::create.");
  }
}

/* This method reads the parameters in */
template <typename Modus>
void Experiment<Modus>::configure(std::list<Parameters>
                                                configuration) {
    assign_params(&configuration);
    warn_wrong_params(&configuration);
    print_startup();
    /* reducing cross section according to number of test particle */
    if (parameters_.testparticles > 1) {
      printf("IC test particle: %i\n", parameters_.testparticles);
      parameters_.cross_section = parameters_.cross_section / parameters_.testparticles;
      printf("Elastic cross section: %g mb\n", parameters_.cross_section);
    }
}

/* This allows command line arguments to override default. */
template <typename Modus>
void Experiment<Modus>::commandline_arg(int steps) {
  steps_ = steps;
}

/* This method reads the particle type and cross section information
 * and does the initialization of the system (fill the particles map)
 */
template <typename Modus>
void Experiment<Modus>::initialize(const char *path) {
    /* Ensure safe allocation */
    delete particles_;
    delete cross_sections_;
    /* Allocate private pointer members */
    particles_ = new Particles;
    cross_sections_ = new CrossSections;
    /* Set the seed_ for the random number generator */
    srand48(seed_);
    /* Read in particle types used in the simulation */
    input_particles(particles_, path);
    /* Read in the particle decay modes */
    input_decaymodes(particles_, path);
    /* Set the default elastic collision cross section */
    cross_sections_->add_elastic_parameter(parameters_.cross_section);
    /* Sample particles according to the initial conditions */
    modus_.initial_conditions(particles_, parameters_);
    /* Save the initial energy in the system for energy conservation checks */
    energy_initial_ = energy_total(particles_);
    /* Print output headers */
    write_measurements_header(*particles_);
    print_header();
    /* Write out the initial momenta and positions of the particles */
    write_particles(*particles_);
}

/* This is the loop over timesteps, carrying out collisions and decays
 * and propagating particles
 */
template <typename Modus>
void Experiment<Modus>::run_time_evolution() {
    modus_.sanity_check(particles_);
    std::list<int> collision_list, decay_list;
    size_t interactions_total = 0, previous_interactions_total = 0,
    interactions_this_interval = 0;
    size_t rejection_conflict = 0;
    int resonances = 0, decays = 0;
    print_measurements(*particles_, interactions_total,
                       interactions_this_interval, energy_initial_,
                       time_start_);
    for (int step = 0; step < steps_; step++) {
        /* Check resonances for decays */
        check_decays(particles_, &decay_list, parameters_.eps);
        /* Do the decays */
        if (!decay_list.empty()) {
            decays += decay_list.size();
            interactions_total = decay_particles(particles_,
                                       &decay_list, interactions_total);
        }
        /* fill collision table by cells */
        modus_.check_collision_geometry(particles_, cross_sections_,
                                 &collision_list, &rejection_conflict, parameters_);
        /* particle interactions */
        if (!collision_list.empty()) {
            printd_list(collision_list);
            interactions_total = collide_particles(particles_, &collision_list,
                                             interactions_total, &resonances);
        }
        modus_.propagate(particles_, parameters_);
        /* physics output during the run */
        if (step > 0 && (step + 1) % output_interval_ == 0) {
            interactions_this_interval = interactions_total
            - previous_interactions_total;
            previous_interactions_total = interactions_total;
            print_measurements(*particles_, interactions_total,
                              interactions_this_interval, energy_initial_,
                              time_start_);
            printd("Resonances: %i Decays: %i\n", resonances, decays);
            printd("Ignored collisions %zu\n", rejection_conflict);
            /* save evolution data */
            write_measurements(*particles_, interactions_total,
                               interactions_this_interval, resonances, decays,
                               rejection_conflict);
            write_vtk(*particles_);
        }
    }
        /* Guard against evolution */
        if (likely(steps_ > 0)) {
            /* if there are no particles no interactions happened */
            if (likely(!particles_->empty())) {
             print_tail(time_start_, interactions_total * 2
                        / particles_->time() / particles_->size());
            } else {
             print_tail(time_start_, 0);
             printf("Total ignored collisions: %zu\n", rejection_conflict);
            }
        }
}

template<typename Modus>
void Experiment<Modus>::assign_params(std::list<Parameters> *configuration) {
    modus_.assign_params(configuration); // TODO: fold into constructor

    bool match = false;
    std::list<Parameters>::iterator i = configuration->begin();
    while (i != configuration->end()) {
        char *key = i->key();
        char *value = i->value();
        printd("Looking for match %s %s\n", key, value);
        /* integer values */
        if (strcmp(key, "STEPS") == 0) {
            steps_ = (abs(atoi(value)));
            match = true;
        }
        if (strcmp(key, "RANDOMSEED") == 0) {
            /* negative seed_ means random startup value */
            if (atol(value) > 0)
                seed_ = (atol(value));
            else
                seed_ = (time(NULL));
            match = true;
        }
        if (strcmp(key, "UPDATE") == 0) {
            output_interval_ = (abs(atoi(value)));
            match = true;
        }
        if (strcmp(key, "TESTPARTICLES") == 0) {
            parameters_.testparticles = (abs(atoi(value)));
            match = true;
        }
        /* double or float values */
        if (strcmp(key, "EPS") == 0) {
            parameters_.eps = (fabs(atof(value)));
            match = true;
        }
        if (strcmp(key, "SIGMA") == 0) {
            parameters_.cross_section = (fabs(atof(value)));
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
template <typename Modus>
void Experiment<Modus>::print_startup() {
    printf("Elastic cross section: %g mb\n", parameters_.cross_section);
    printf("Using temporal stepsize: %g fm/c\n", parameters_.eps);
    printf("Maximum number of steps: %i \n", steps_);
    printf("Random number seed: %" PRId64 "\n", seed_);
    modus_.print_startup();
}

/* calculates the total energy in the system from zero component of
 * all momenta of particles
 * XXX should be expanded to all quantum numbers of interest */
template <typename Modus>
float Experiment<Modus>::energy_total(Particles *particles) {
    float energy_sum = 0.0;
    for (auto i = particles->begin(); i != particles->end(); ++i) {
         energy_sum += i->second.momentum().x0();
     }
    return energy_sum;
}

/* set the timer to the actual time in nanoseconds precision */
template <typename Modus>
timespec inline Experiment<Modus>::set_timer_start(void) {
    timespec time;
    clock_gettime(&time);
    return time;
}

/* Tear down everything */
template <typename Modus>
void Experiment<Modus>::end() {
    delete particles_;
    delete cross_sections_;
}

template <typename Modus>
void Experiment<Modus>::run(std::string path) {
  for (int j = 1; j < nevents_; j++) {
    write_oscar_header();
    initialize(path.c_str());
    /* the time evolution of the relevant subsystem */
    run_time_evolution();
  }
  end();
}
