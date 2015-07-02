/*
 *
 *    Copyright (c) 2012-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/experiment.h"

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
#include "include/cxx14compat.h"
#include "include/decayactionsfinder.h"
#include "include/decayactionsfinderdilepton.h"
#include "include/density.h"
#include "include/forwarddeclarations.h"
#include "include/grid.h"
#include "include/listmodus.h"
#include "include/logging.h"
#include "include/macros.h"
#include "include/pauliblocking.h"
#include "include/potentials.h"
#include "include/propagation.h"
#include "include/random.h"
#include "include/scatteractionsfinder.h"
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
   * \li \key Sphere for calculations of the expansion of a thermalized sphere.
   * See
   *     \ref \SphereModus
   * \li \key Box for infinite matter calculation in a rectangular box. See \ref
   *     \BoxModus
   * \li \key List for given external particle list. See \ref
   *     \ListModus
   */

  /*!\Userguide
   * \page input_modi_ Modi
   * \li \subpage input_modi_collider_
   * \li \subpage input_modi_sphere_
   * \li \subpage input_modi_box_
   * \li \subpage input_modi_list_
   */
  const std::string modus_chooser = config.take({"General", "Modus"});
  log.info() << "Modus for this calculation: " << modus_chooser;

  // remove config maps of unused Modi
  config["Modi"].remove_all_but(modus_chooser);

  typedef std::unique_ptr<ExperimentBase> ExperimentPointer;
  if (modus_chooser.compare("Box") == 0) {
    if (config.has_value({"Potentials"})) {
      log.error() << "Box modus does not work with potentials for now: "
                  << "periodic boundaries are not taken into account "
                  << "in the density calculation";
    }
    return ExperimentPointer(new Experiment<BoxModus>(config));
  } else if (modus_chooser.compare("List") == 0) {
    return ExperimentPointer(new Experiment<ListModus>(config));
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
 * \key Gaussian_Sigma (float, optional, default 1.0): \n
 * Width [fm] of gaussians that represent Wigner density of particles.
 *
 * \key Gaussian_Cutoff (float, optional, default 4.0)
 * Distance in sigma at which gaussian is considered 0.
 *
 * \page input_output_options_ Output
 * \key Output_Interval (float, required): \n
 * Defines the period of intermediate output of the status of the simulated
 * system in Standard Output and other output formats which support this
 * functionality.
 *
 * \key Gaussian_Sigma (float, optional, default 1.0): \n
 * Width [fm] of gaussians that represent Wigner density of particles.
 *
 * \key Gaussian_Cutoff (float, optional, default 4.0)
 * Distance in sigma at which gaussian is considered 0.
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

  // The clock initializers are only read here and taken later when
  // assigning initial_clock_.
  return {{0.0f, config.read({"General", "Delta_Time"})},
          config.take({"Output", "Output_Interval"}),
          config.take({"General", "Testparticles"}, 1),
          config.take({"General", "Gaussian_Sigma"}, 1.0),
          config.take({"General", "Gauss_Cutoff_In_Sigma"}, 4.0)};
}
}  // unnamed namespace

/**
 * Creates a verbose textual description of the setup of the Experiment.
 */
template <typename Modus>
std::ostream &operator<<(std::ostream &out, const Experiment<Modus> &e) {
  out << "Starting with temporal stepsize: "
      << e.parameters_.timestep_duration() << " fm/c\n";
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
 *
 * \key Density_Type (int, optional, default = 0): \n
 * 0 - net baryon density
 * 1 - baryonic isospin density
 *
 * \key Force_Decays_At_End (bool, optional, default = true): \n
 * true - force all resonances to decay after last timestep \n
 * false - don't force decays (final output can contain resonances)
 */
template <typename Modus>
Experiment<Modus>::Experiment(Configuration config)
    : parameters_(create_experiment_parameters(config)),
      modus_(config["Modi"], parameters_),
      particles_(),
      nevents_(config.take({"General", "Nevents"})),
      end_time_(config.take({"General", "End_Time"})),
      delta_time_startup_(config.take({"General", "Delta_Time"})),
      force_decays_(
          config.take({"Collision_Term", "Force_Decays_At_End"}, true)) {
  const auto &log = logger<LogArea::Experiment>();
  log.info() << *this;

  if (config.take({"Collision_Term", "Decays"}, true)) {
    action_finders_.emplace_back(new DecayActionsFinder());
  }
  if (config.take({"Collision_Term", "Dileptons"}, false)) {
  action_finders_.emplace_back(new DecayActionsFinderDilepton(bf::path(bf::current_path()),config.take({"Collision_Term", "Dileptons"}, "dileptons"));
} // worst hack ever for config value
  if (config.take({"Collision_Term", "Collisions"}, true)) {
    action_finders_.emplace_back(new ScatterActionsFinder(config, parameters_));
  }
  if (config.has_value({"Collision_Term", "Pauli_Blocking"})) {
    log.info() << "Pauli blocking is ON.";
    pauli_blocker_ = make_unique<PauliBlocker>(
        config["Collision_Term"]["Pauli_Blocking"], parameters_);
  }

  if (config.has_value({"Potentials"})) {
    log.info() << "Potentials are ON.";
    // potentials need testparticles and gaussian sigma from parameters_
    potentials_ = make_unique<Potentials>(config["Potentials"], parameters_);
  }

  dens_type_ =
      config.take({"Output", "Density", "Density_Type"}, DensityType::particle);
  log.info() << "Density type written to headers: " << dens_type_;

  // Create lattices
  if (config.has_value({"Lattice"})) {
    // Take lattice properties from config to assign them to all lattices
    const std::array<float, 3> l = config.take({"Lattice", "Sizes"});
    const std::array<int, 3> n = config.take({"Lattice", "CellNumber"});
    const std::array<float, 3> origin = config.take({"Lattice", "Origin"});
    const bool periodic = config.take({"Lattice", "Periodic"});
    dens_type_lattice_printout_ = static_cast<DensityType>(
         config.take({"Lattice", "Printout", "Density"},
                     static_cast<int>(DensityType::none)));
    /* Create baryon and isospin density lattices regardless of config
       if potentials are on. This is because they allow to compute
       potentials faster */
    if (potentials_) {
      jmu_B_lat_ = make_unique<DensityLattice>(l, n, origin, periodic,
                                            LatticeUpdate::EveryTimestep);
      jmu_I3_lat_ = make_unique<DensityLattice>(l, n, origin, periodic,
                                              LatticeUpdate::EveryTimestep);
    } else {
      if (dens_type_lattice_printout_ == DensityType::baryon) {
        jmu_B_lat_ = make_unique<DensityLattice>(l, n, origin, periodic,
                                                  LatticeUpdate::AtOutput);
      }
      if (dens_type_lattice_printout_ == DensityType::baryonic_isospin) {
        jmu_I3_lat_ = make_unique<DensityLattice>(l, n, origin, periodic,
                                             LatticeUpdate::AtOutput);
      }
    }
    if (dens_type_lattice_printout_ != DensityType::none &&
        dens_type_lattice_printout_ != DensityType::baryonic_isospin &&
        dens_type_lattice_printout_ != DensityType::baryon) {
        jmu_custom_lat_ = make_unique<DensityLattice>(l, n, origin,
                                          periodic, LatticeUpdate::AtOutput);
    }
  }
}

const std::string hline(80, '-');

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
  conserved_initial_ = QuantumNumbers(particles_);
  /* Print output headers */
  log.info() << hline;
  log.info() << " Time       <Ediff>      <pdiff>  <scattrate>    <scatt>  "
                "<particles>   <timing>";
  log.info() << hline;
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

  std::ostringstream ss;
  ss << field<5> << time
     << field<12, 3> << difference.momentum().x0()
     << field<12, 3> << difference.momentum().abs3()
     << field<12, 3> << (scatterings_total
                          ? scatterings_total * 2 / (particles.size() * time)
                          : 0.)
     << field<10, 3> << scatterings_this_interval
     << field<12, 3> << particles.size()
     << field<10, 3> << elapsed_seconds;
  return ss.str();
}

template <typename Modus>
void Experiment<Modus>::perform_actions(ActionList &actions,
                                        size_t &interactions_total,
                                        size_t &total_pauli_blocked) {
  const auto &log = logger<LogArea::Experiment>();
  if (!actions.empty()) {
    const auto particles_before_actions = particles_.copy_to_vector();
    for (const auto &action : actions) {
      if (action->is_valid(particles_)) {
        const ParticleList incoming_particles = action->incoming_particles();
        action->generate_final_state();
        ProcessType process_type = action->get_type();
        log.debug("Process Type is: ", process_type);
        if (pauli_blocker_ &&
            action->is_pauli_blocked(particles_, *pauli_blocker_.get())) {
          total_pauli_blocked++;
          continue;
        }
        action->perform(&particles_, interactions_total);
        const ParticleList outgoing_particles = action->outgoing_particles();
        // Calculate Eckart rest frame density at the interaction point
        const FourVector r_interaction = action->get_interaction_point();
        constexpr bool compute_grad = false;
        const double rho =
            rho_eckart(r_interaction.threevec(), particles_before_actions,
                       parameters_, dens_type_, compute_grad).first;
        /*!\Userguide
         * \page collisions_output_in_box_modus_ Collision output in box modus
         * \note When SMASH is running in the box modus, particle coordinates
         * in the collision output can be out of the box. This is not an error.
         * Box boundary conditions are intentionally not imposed before
         * collision output to allow unambiguous finding of the interaction
         * point.
         * <I>Example</I>: two particles in the box have x coordinates 0.1 and
         * 9.9 fm, while box L = 10 fm. Suppose these particles collide.
         * For calculating collision the first one is wrapped to 10.1 fm.
         * Then output contains coordinates of 9.9 fm and 10.1 fm.
         * From this one can infer interaction point at x = 10 fm.
         * Were boundary conditions imposed before output,
         * their x coordinates would be 0.1 and 9.9 fm and interaction point
         * position could be either at 10 fm or at 5 fm.
         */
        for (const auto &output : outputs_) {
          output->at_interaction(incoming_particles, outgoing_particles, rho,
                                 action->raw_weight_value(), process_type);
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
}

/* This is the loop over timesteps, carrying out collisions and decays
 * and propagating particles. */
template <typename Modus>
void Experiment<Modus>::run_time_evolution(const int evt_num) {
  const auto &log = logger<LogArea::Experiment>();
  modus_.impose_boundary_conditions(&particles_);
  size_t interactions_total = 0, previous_interactions_total = 0,
         interactions_this_interval = 0, total_pauli_blocked = 0;
  log.info() << format_measurements(
      particles_, interactions_total, interactions_this_interval,
      conserved_initial_, time_start_, parameters_.labclock.current_time());

  while (!(++parameters_.labclock > end_time_)) {
    // vector is likely the best container type here. Because std::sort requires
    // random access iterators. Any linked data structure (e.g. list) thus
    // requires a less efficient sort algorithm.
    std::vector<ActionPtr> actions;

    /* (1.a) Create grid. */
    const auto &grid =
        // TODO(mkretz): avoid the copy. Grid could construct from Particles
        // directly.
        modus_.create_grid(particles_, parameters_.testparticles);
    /* (1.b) Iterate over cells and find actions. */
    grid.iterate_cells([&](const ParticleList &search_list) {
                         for (const auto &finder : action_finders_) {
                           actions += finder->find_possible_actions(
                               search_list, parameters_.timestep_duration());
                         }
                       },
                       [&](const ParticleList &search_list,
                           const ParticleList &neighbors_list) {
                         for (const auto &finder : action_finders_) {
                           actions += finder->find_possible_actions(
                               search_list, neighbors_list,
                               parameters_.timestep_duration());
                         }
                       });
    /* (1.c) Sort action list by time. */
    std::sort(actions.begin(), actions.end(),
              [](const ActionPtr &a, const ActionPtr &b) { return *a < *b; });

    /* (2) Perform actions. */
    perform_actions(actions, interactions_total, total_pauli_blocked);
    modus_.impose_boundary_conditions(&particles_);

    /* (3) Do propagation. */
    if (potentials_) {
      update_density_lattice(jmu_B_lat_.get(), LatticeUpdate::EveryTimestep,
                       DensityType::baryon, parameters_, particles_);
      update_density_lattice(jmu_I3_lat_.get(), LatticeUpdate::EveryTimestep,
                       DensityType::baryonic_isospin, parameters_, particles_);
      propagate(&particles_, parameters_, *potentials_);
    } else {
      propagate_straight_line(&particles_, parameters_);
    }
    modus_.impose_boundary_conditions(&particles_, outputs_);

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
      // Update lattices for output
      const LatticeUpdate lat_upd = LatticeUpdate::AtOutput;
      update_density_lattice(jmu_B_lat_.get(), lat_upd, DensityType::baryon,
                             parameters_, particles_);
      update_density_lattice(jmu_I3_lat_.get(), lat_upd,
                DensityType::baryonic_isospin, parameters_, particles_);
      update_density_lattice(jmu_custom_lat_.get(), lat_upd,
                dens_type_lattice_printout_, parameters_, particles_);
      /* save evolution data */
      for (const auto &output : outputs_) {
        output->at_intermediate_time(particles_, evt_num, parameters_.labclock);
        output->thermodynamics_output(particles_, parameters_);
        switch (dens_type_lattice_printout_) {
          case DensityType::baryon:
            output->thermodynamics_output(std::string("rhoB"), *jmu_B_lat_,
                                                                     evt_num);
            break;
          case DensityType::baryonic_isospin:
            output->thermodynamics_output(std::string("rhoI3"), *jmu_I3_lat_,
                                                                     evt_num);
            break;
          case DensityType::none:
            break;
          default:
            output->thermodynamics_output(std::string("rho"), *jmu_custom_lat_,
                                                                      evt_num);
        }
      }
    }
    // Check conservation of conserved quantities if potentials are off.
    // If potentials are on then momentum is conserved only in average
    if (!potentials_) {
      std::string err_msg = conserved_initial_.report_deviations(particles_);
      if (!err_msg.empty()) {
        log.error() << err_msg;
        throw std::runtime_error("Violation of conserved quantities!");
      }
    }
  }

  if (pauli_blocker_) {
    log.info("Collisions: pauliblocked/total = ", total_pauli_blocked, "/",
             interactions_total);
  }

  if (force_decays_) {
    // at end of time evolution: force all resonances to decay
    size_t interactions_old;
    do {
      std::vector<ActionPtr> actions;
      interactions_old = interactions_total;
      /* Find actions. */
      for (const auto &finder : action_finders_) {
        actions += finder->find_final_actions(particles_);
      }
      /* Perform actions. */
      perform_actions(actions, interactions_total, total_pauli_blocked);
      // loop until no more decays occur
    } while (interactions_total > interactions_old);

    /* Do one final propagation step. */
    if (potentials_) {
      propagate(&particles_, parameters_, *potentials_);
    } else {
      propagate_straight_line(&particles_, parameters_);
    }
    modus_.impose_boundary_conditions(&particles_, outputs_);
  }

  // make sure the experiment actually ran (note: we should compare this
  // to the start time, but we don't know that. Therefore, we check that
  // the time is positive, which should heuristically be the same).
  if (likely(parameters_.labclock > 0)) {
    log.info() << hline;
    log.info() << "Time real: " << SystemClock::now() - time_start_;
    /* if there are no particles no interactions happened */
    log.info() << "Final scattering rate: "
               << (particles_.is_empty() ? 0 : (interactions_total * 2 /
                                                particles_.time() /
                                                particles_.size()))
               << " [fm-1]";
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
