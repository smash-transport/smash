/*
 *
 *    Copyright (c) 2012-2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */


#include <cinttypes>
#include <cstdlib>
#include <ctime>
#include <list>
#include <string>
#include <algorithm>

#include "include/boxmodus.h"
#include "include/collidermodus.h"
#include "include/configuration.h"
#include "include/experiment.h"
#include "include/macros.h"
#include "include/nucleusmodus.h"
#include "include/outputroutines.h"
#include "include/particlesoutput.h"
#include "include/random.h"
#include "include/time.h"
#include "include/vtkoutput.h"

#include <boost/filesystem.hpp>

/* #include "include/spheremodus.h" */

namespace Smash {

namespace bf = boost::filesystem;

/* ExperimentBase carries everything that is needed for the evolution */
std::unique_ptr<ExperimentBase> ExperimentBase::create(Configuration &config) {
  const std::string modus_chooser = config.take({"General", "MODUS"});
  printf("Modus for this calculation: %s\n", modus_chooser.c_str());

  // remove config maps of unused Modi
  config["Modi"].remove_all_but(modus_chooser);

  typedef std::unique_ptr<ExperimentBase> ExperimentPointer;
  if (modus_chooser.compare("Box") == 0) {
    return ExperimentPointer(new Experiment<BoxModus>(config));
  } else if (modus_chooser.compare("Collider") == 0) {
    return ExperimentPointer(new Experiment<ColliderModus>(config));
  } else if (modus_chooser.compare("Nucleus") == 0) {
    return ExperimentPointer(new Experiment<NucleusModus>(config));
  } else {
    throw InvalidModusRequest("Invalid Modus (" + modus_chooser +
                              ") requested from ExperimentBase::create.");
  }
}

namespace {
ExperimentParameters create_experiment_parameters(Configuration &config) {
  const int testparticles = config.take({"General", "TESTPARTICLES"});
  float cross_section = config.take({"General", "SIGMA"});

  /* reducing cross section according to number of test particle */
  if (testparticles > 1) {
    printf("IC test particle: %i\n", testparticles);
    cross_section /= testparticles;
    printf("Elastic cross section: %g mb\n", cross_section);
  }

  return {config.take({"General", "EPS"}), cross_section, testparticles};
}
}  // unnamed namespace

template <typename Modus>
Experiment<Modus>::Experiment(Configuration &config)
    : parameters_(create_experiment_parameters(config)),
      modus_(config["Modi"], parameters_),
      particles_{config.take({"particles"}), config.take({"decaymodes"})},
      cross_sections_(parameters_.cross_section),
      nevents_(config.take({"General", "NEVENTS"})),
      steps_(config.take({"General", "STEPS"})),
      output_interval_(config.take({"General", "UPDATE"})) {
  int64_t seed_ = config.take({"General", "RANDOMSEED"});
  if (seed_ < 0) {
    seed_ = time(nullptr);
  }
  Random::set_seed(seed_);

  print_startup(seed_);
}

/* This method reads the particle type and cross section information
 * and does the initialization of the system (fill the particles map)
 */
template <typename Modus>
void Experiment<Modus>::initialize(const bf::path &/*path*/) {
  cross_sections_.reset();
  particles_.reset();

  /* Sample particles according to the initial conditions */
  modus_.initial_conditions(&particles_, parameters_);
  /* Save the initial energy in the system for energy conservation checks */
  energy_initial_ = energy_total(&particles_);
  /* Print output headers */
  print_header();
  /* Write out the initial momenta and positions of the particles */
  for (auto &output : outputs_) {
    output->write_state(particles_);
  }
}


/* This is the loop over timesteps, carrying out collisions and decays
 * and propagating particles. */
template <typename Modus>
void Experiment<Modus>::run_time_evolution() {
  modus_.sanity_check(&particles_);
  size_t interactions_total = 0, previous_interactions_total = 0,
         interactions_this_interval = 0;
  print_measurements(particles_, interactions_total,
                     interactions_this_interval, energy_initial_, time_start_);

  for (int step = 0; step < steps_; step++) {
    std::vector<ActionPtr> actions;

    /* (1.a) Find possible decays. */
    actions += decay_finder_.find_possible_actions(&particles_, parameters_);
    /* (1.b) Find possible collisions. */
    actions += scatter_finder_.find_possible_actions(&particles_, parameters_,
                                                     &cross_sections_);
    /* (1.c) Sort action list by time. */
    std::sort (actions.begin(), actions.end());

    /* (2.a) Perform actions. */
    if (!actions.empty()) {
      for (auto act = actions.begin(); act != actions.end(); ++act) {
        (*act)->perform (&particles_, interactions_total);
      }
      actions.clear();
      printd("Action list done.\n");
    }

    /* (3) Do propagation. */
    modus_.propagate(&particles_, parameters_);

    /* (4) Physics output during the run. */
    if (step > 0 && (step + 1) % output_interval_ == 0) {
      interactions_this_interval =
          interactions_total - previous_interactions_total;
      previous_interactions_total = interactions_total;
      print_measurements(particles_, interactions_total,
                         interactions_this_interval, energy_initial_,
                         time_start_);
      /* save evolution data */
      for (auto &output : outputs_) {
        output->write_state(particles_);
      }
    }
  }
  /* Guard against evolution */
  if (likely(steps_ > 0)) {
    /* if there are no particles no interactions happened */
    if (likely(!particles_.empty())) {
      print_tail(time_start_, interactions_total * 2 / particles_.time() /
                                  particles_.size());
    } else {
      print_tail(time_start_, 0);
    }
  }
}

/* print_startup - console output on startup of general parameters */
template <typename Modus>
void Experiment<Modus>::print_startup(int64_t seed) {
  printf("Elastic cross section: %g mb\n", parameters_.cross_section);
  printf("Using temporal stepsize: %g fm/c\n", parameters_.eps);
  printf("Maximum number of steps: %i \n", steps_);
  printf("Random number seed: %" PRId64 "\n", seed);
  modus_.print_startup();
}

/* calculates the total energy in the system from zero component of
 * all momenta of particles
 * XXX should be expanded to all quantum numbers of interest */
template <typename Modus>
float Experiment<Modus>::energy_total(Particles *particles) {
  float energy_sum = 0.0;
  for (const ParticleData &data : particles->data()) {
    energy_sum += data.momentum().x0();
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

template <typename Modus>
void Experiment<Modus>::run(const bf::path &path) {
  outputs_.emplace_back(new ParticlesOutput(path));
  outputs_.emplace_back(new VtkOutput(path));

  /* Write the header of OSCAR data output file */
  write_oscar_header();
  for (int j = 0; j < nevents_; j++) {
    initialize(path);
    /* Write the initial data block of the event */
    write_oscar_event_block(&particles_, 0, particles_.size(), j + 1);
    /* the time evolution of the relevant subsystem */
    run_time_evolution();
    /* Write the final data block of the event */
    write_oscar_event_block(&particles_, particles_.size(), 0, j + 1);
  }
}

}  // namespace Smash
