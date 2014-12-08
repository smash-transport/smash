/*
 *
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */


#include <algorithm>
#include <cinttypes>
#include <cstdlib>
#include <list>
#include <string>
#include <vector>

#include "include/action.h"
#include "include/boxmodus.h"
#include "include/clock.h"
#include "include/collidermodus.h"
#include "include/configuration.h"
#include "include/experiment.h"
#include "include/forwarddeclarations.h"
#include "include/grid.h"
#include "include/logging.h"
#include "include/macros.h"
#include "include/random.h"
#include "include/spheremodus.h"

namespace std {
/**
 * Print time span in a human readable way:
 * time < 10 min => seconds
 * 10 min < time < 3 h => minutes
 * time > 3h => hours
 *
 * \note This operator has to be in the \c std namespace for argument dependent
 * lookup to find it. If it were in the Smash namespace then the code would not
 * compile since none of its arguments is a type from the Smash namespace.
 */
template <typename T, typename Ratio>
static ostream &operator<<(ostream &out,
                           const chrono::duration<T, Ratio> &seconds) {
  using namespace chrono;
  using Seconds = duration<float>;
  using Minutes = duration<float, std::ratio<60>>;
  using Hours = duration<float, std::ratio<60 * 60>>;
  constexpr Minutes threshold_for_minutes{10};
  constexpr Hours threshold_for_hours{3};
  if (seconds < threshold_for_minutes) {
    return out << Seconds(seconds).count() << " [s]";
  }
  if (seconds < threshold_for_hours) {
    return out << Minutes(seconds).count() << " [min]";
  }
  return out << Hours(seconds).count() << " [h]";
}
}  // namespace std

namespace Smash {

/* ExperimentBase carries everything that is needed for the evolution */
std::unique_ptr<ExperimentBase> ExperimentBase::create(Configuration config) {
  const auto &log = logger<LogArea::Experiment>();
  log.trace() << source_location;
  /*!\Userguide
   * \page input_general_ General
   * \key Modus (string, required): \n
   * Selects a modus for the calculation, e.g.\ infinite matter
   * calculation, collision of two particles or collision of nuclei. The modus
   * will be configured in \ref input_modi_. Recognized values are:
   * \li \key Collider for collisions of nuclei or compound objects. See \ref
   *     \ColliderModus
   * \li \key Sphere for calculations of the expansion of a thermalized sphere. See
   *     \ref \SphereModus
   * \li \key Box for infinite matter calculation in a rectangular box. See \ref
   *     \BoxModus
   */

  /*!\Userguide
   * \page input_modi_ Modi
   * \li \subpage input_modi_collider_
   * \li \subpage input_modi_sphere_
   * \li \subpage input_modi_box_
   */
  const std::string modus_chooser = config.take({"General", "Modus"});
  log.info() << "Modus for this calculation: " << modus_chooser;

  // remove config maps of unused Modi
  config["Modi"].remove_all_but(modus_chooser);

  typedef std::unique_ptr<ExperimentBase> ExperimentPointer;
  if (modus_chooser.compare("Box") == 0) {
    return ExperimentPointer(new Experiment<BoxModus>(config));
  } else if (modus_chooser.compare("Collider") == 0) {
      return ExperimentPointer(new Experiment<ColliderModus>(config));
  } else if (modus_chooser.compare("Sphere") == 0) {
      return ExperimentPointer(new Experiment<SphereModus>(config));
  } else {
    throw InvalidModusRequest("Invalid Modus (" + modus_chooser +
                              ") requested from ExperimentBase::create.");
  }
}

namespace {
/*!\Userguide
 * \page input_general_ General
 * \key Delta_Time (float, required): \n
 * Time step for the calculation, in fm/c.
 *
 * \key Testparticles (int, optional, default = 1): \n
 * How many test particles per real particles should be simulated.
 *
 * \key Gaussian_Sigma (float, required): \n
 * Width of gaussians that represent Wigner density of particles.
 *
 * \page input_output_options_ Output
 * \key Output_Interval (float, required): \n
 * Defines the period of intermediate output of the status of the simulated
 * system in Standard Output and other output formats which support this
 * functionality.
 */
/** Gathers all general Experiment parameters
 *
 * \param[in, out] config Configuration element
 * \return The ExperimentParameters struct filled with values from the
 * Configuration
 */
ExperimentParameters create_experiment_parameters(Configuration config) {
  const auto &log = logger<LogArea::Experiment>();
  log.trace() << source_location;

  int testparticles;
  if (config.has_value({"General", "Testparticles"})) {
    testparticles = config.take({"General", "Testparticles"});
  } else {
    testparticles = 1;
  }
  // The clock initializers are only read here and taken later when
  // assigning initial_clock_.
  return {{0.0f, config.read({"General", "Delta_Time"})},
           config.take({"Output", "Output_Interval"}),
           testparticles, config.take({"General", "Gaussian_Sigma"})};
  }
}  // unnamed namespace

/**
 * Creates a verbose textual description of the setup of the Experiment.
 */
template <typename Modus>
std::ostream &operator<<(std::ostream &out, const Experiment<Modus> &e) {
  out << "Starting with temporal stepsize: " << e.parameters_.timestep_duration()
      << " fm/c\n";
  out << "End time: " << e.end_time_ << " fm/c\n";
  out << e.modus_;
  return out;
}

/*!\Userguide
 * \page input_general_
 * \key End_Time (float, required): \n
 * The time after which the evolution is stopped. Note
 * that the starting time depends on the chosen Modus.
 *
 * \key Randomseed (int64_t, required): \n
 * Initial seed for the random number generator. If this is
 * negative, the program starting time is used.
 *
 * \key Nevents (int, required): \n
 * Number of events to calculate.
 *
 * \page input_collision_term_ Collision_Term
 * \key Decays (bool, optional, default = true): \n
 * true - decays are enabled\n
 * false - disable all decays
 *
 * \key Collisions (bool, optional, default = true): \n
 * true - collisions are enabled\n
 * false - all collisions are disabled
 */
template <typename Modus>
Experiment<Modus>::Experiment(Configuration config)
    : parameters_(create_experiment_parameters(config)),
      modus_(config["Modi"], parameters_),
      particles_(),
      nevents_(config.take({"General", "Nevents"})),
      end_time_(config.take({"General", "End_Time"})),
      delta_time_startup_(config.take({"General", "Delta_Time"})) {
  const auto &log = logger<LogArea::Experiment>();
  int64_t seed_ = config.take({"General",  "Randomseed"});
  if (seed_ < 0) {
    // Seed with a real random value, if available
    std::random_device rd;
    seed_ = rd();
  }
  Random::set_seed(seed_);
  log.info() << "Random number seed: " << seed_;
  log.info() << *this;

  if (!config.has_value({"Collision_Term", "Decays"})
      || config.take({"Collision_Term", "Decays"})) {
    action_finders_.emplace_back(new DecayActionsFinder(parameters_));
  }
  if (!config.has_value({"Collision_Term", "Collisions"})
      || config.take({"Collision_Term", "Collisions"})) {
    action_finders_.emplace_back(new ScatterActionsFinder(config, parameters_));
  }
}

/* This method reads the particle type and cross section information
 * and does the initialization of the system (fill the particles map)
 */
template <typename Modus>
void Experiment<Modus>::initialize_new_event() {
  const auto &log = logger<LogArea::Experiment>();
  particles_.reset();

  /* Sample particles according to the initial conditions */
  float start_time = modus_.initial_conditions(&particles_, parameters_);

  // reset the clock:
  Clock clock_for_this_event(start_time, delta_time_startup_);
  parameters_.labclock = std::move(clock_for_this_event);

  /* Save the initial conserved quantum numbers and total momentum in
   * the system for conservation checks */
  conserved_initial_.count_conserved_values(particles_);
  /* Print output headers */
  log.info() << "--------------------------------------------------------------------------------";
  log.info() << " Time       <Ediff>      <pdiff>  <scattrate>    <scatt>  <particles>   <timing>";
  log.info() << "--------------------------------------------------------------------------------";
}

static std::string format_measurements(const Particles &particles,
                                       size_t scatterings_total,
                                       size_t scatterings_this_interval,
                                       const QuantumNumbers &conserved_initial,
                                       SystemTimePoint time_start,
                                       double time) {
  SystemTimeSpan elapsed_seconds = SystemClock::now() - time_start;

  QuantumNumbers current_values(particles);
  QuantumNumbers difference = conserved_initial - current_values;

  char buffer[81];
  if (likely(time > 0))
    snprintf(buffer, 81, "%6.2f %12g %12g %12g %10zu %12zu %12g", time,
             difference.momentum().x0(), difference.momentum().abs3(),
             scatterings_total * 2 / (particles.size() * time),
             scatterings_this_interval, particles.size(),
             elapsed_seconds.count());
  else
    snprintf(buffer, 81, "%+6.2f %12g %12g %12g %10i %12zu %12g", time,
             difference.momentum().x0(), difference.momentum().abs3(), 0.0, 0,
             particles.size(), elapsed_seconds.count());
  return std::string(buffer);
}

/* This is the loop over timesteps, carrying out collisions and decays
 * and propagating particles. */
template <typename Modus>
void Experiment<Modus>::run_time_evolution(const int evt_num) {
  const auto &log = logger<LogArea::Experiment>();
  modus_.sanity_check(&particles_);
  size_t interactions_total = 0, previous_interactions_total = 0,
         interactions_this_interval = 0;
  log.info() << format_measurements(
      particles_, interactions_total, interactions_this_interval,
      conserved_initial_, time_start_, parameters_.labclock.current_time());

  while (!(++parameters_.labclock > end_time_)) {
    std::vector<ActionPtr> actions;  // TODO: a std::list might be better suited
                                     // for the task: lots of appending, then
                                     // sorting and finally a single linear
                                     // iteration

    const auto &grid = modus_.create_grid(
        ParticleList{particles_.data().begin(), particles_.data().end()});
    grid.iterate_cells([&](
        const ParticleList &search_list,  // a list of particles where each pair
                                          // needs to be tested for possible
                                          // interaction
        const std::vector<const ParticleList *> &
            neighbors_list  // a list of particles that need to be tested
                            // against particles in search_list for possible
                            // interaction
        ) {
      for (const auto &finder : action_finders_) {
        actions += finder->find_possible_actions(search_list, neighbors_list);
      }
    });

    /* (1.c) Sort action list by time. */
    std::sort(actions.begin(), actions.end(),
              [](const ActionPtr &a, const ActionPtr &b) { return *a < *b; });

    /* (2.a) Perform actions. */
    if (!actions.empty()) {
      for (const auto &action : actions) {
        if (action->is_valid(particles_)) {
          const ParticleList incoming_particles = action->incoming_particles();
          action->perform(&particles_, interactions_total);
          const ParticleList outgoing_particles = action->outgoing_particles();
          for (const auto &output : outputs_) {
            output->at_interaction(incoming_particles, outgoing_particles);
          }
          log.debug(~einhard::Green(), "✔ ", action);
        } else {
          log.debug(~einhard::DRed(), "✘ ", action, " (discarded: invalid)");
        }
      }
      actions.clear();
      log.debug(~einhard::Blue(), particles_);
    } else {
      log.debug("no actions performed");
    }

    modus_.sanity_check(&particles_);

    /* (3) Do propagation. */
    modus_.propagate(&particles_, parameters_, outputs_);

    /* (4) Physics output during the run. */
    // if the timestep of labclock is different in the next tick than
    // in the current one, I assume it has been changed already. In that
    // case, I know what the next tick is and I can check whether the
    // output time is crossed within the next tick.
    if (parameters_.need_intermediate_output()) {
      interactions_this_interval =
          interactions_total - previous_interactions_total;
      previous_interactions_total = interactions_total;
      log.info() << format_measurements(
          particles_, interactions_total, interactions_this_interval,
          conserved_initial_, time_start_, parameters_.labclock.current_time());
      /* save evolution data */
      for (const auto &output : outputs_) {
        output->at_intermediate_time(particles_, evt_num, parameters_.labclock);
      }
    }
    // check conservation of conserved quantities:
    std::string err_msg = conserved_initial_.report_deviations(particles_);
    if (!err_msg.empty()) {
      log.error() << err_msg;
      throw std::runtime_error("Violation of conserved quantities!");
    }
  }
  // make sure the experiment actually ran (note: we should compare this
  // to the start time, but we don't know that. Therefore, we check that
  // the time is positive, which should heuristically be the same).
  if (likely(parameters_.labclock > 0)) {
    log.info() << "--------------------------------------------------------------------------------";
    log.info() << "Time real: " << SystemClock::now() - time_start_;
    /* if there are no particles no interactions happened */
    log.info() << "Final scattering rate: "
               << (particles_.empty() ? 0 : (interactions_total * 2 /
                                             particles_.time() /
                                             particles_.size())) << " [fm-1]";
  }
}

template <typename Modus>
void Experiment<Modus>::run() {
  const auto &mainlog = logger<LogArea::Main>();
  for (int j = 0; j < nevents_; j++) {
    mainlog.info() << "Event " << j;

    /* Sample initial particles, start clock, some printout and book-keeping */
    initialize_new_event();

    /* Output at event start */
    for (const auto &output : outputs_) {
      output->at_eventstart(particles_, j);
    }

    /* the time evolution of the relevant subsystem */
    run_time_evolution(j);

    /* Output at event end */
    for (const auto &output : outputs_) {
      output->at_eventend(particles_, j);
    }
  }
}

}  // namespace Smash
