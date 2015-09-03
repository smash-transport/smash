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
#include "include/actions.h"
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
/* Outputs */
#include "include/binaryoutputcollisions.h"
#include "include/binaryoutputparticles.h"
#include "include/densityoutput.h"
#include "include/oscaroutput.h"
#ifdef SMASH_USE_ROOT
#  include "include/rootoutput.h"
#endif
#include "include/vtkoutput.h"

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
  using Seconds = chrono::duration<float>;
  using Minutes = chrono::duration<float, std::ratio<60>>;
  using Hours = chrono::duration<float, std::ratio<60 * 60>>;
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
std::unique_ptr<ExperimentBase> ExperimentBase::create(Configuration config,
                                                       bf::path output_path) {
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
    if (config.has_value({"General", "Time_Step_Mode"}) &&
        config.read({"General", "Time_Step_Mode"}) == TimeStepMode::None) {
      log.error() << "Box modus does not work correctly without time steps for "
                  << "now: periodic boundaries are not taken into account when "
                  << "looking for interactions.";
      throw std::invalid_argument("Can't use box modus without time steps!");
    }
    return ExperimentPointer(new Experiment<BoxModus>(config, output_path));
  } else if (modus_chooser.compare("List") == 0) {
    return ExperimentPointer(new Experiment<ListModus>(config, output_path));
  } else if (modus_chooser.compare("Collider") == 0) {
    return ExperimentPointer(new Experiment<ColliderModus>(config,
                                                           output_path));
  } else if (modus_chooser.compare("Sphere") == 0) {
    return ExperimentPointer(new Experiment<SphereModus>(config, output_path));
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
 * \key Gauss_Cutoff_In_Sigma (float, optional, default 4.0)
 * Distance in sigma at which gaussian is considered 0.
 *
 * \page input_output_options_ Output
 * \key Output_Interval (float, required): \n
 * Defines the period of intermediate output of the status of the simulated
 * system in Standard Output and other output formats which support this
 * functionality.
 *
 * \key Density_Type (string, optional, default = "none"): \n
 * Determines which kind of density is written into the collision files.
 * Possible values:\n
 * \li "hadron"           - total hadronic density \n
 * \li "baryon"           - net baryon density \n
 * \li "baryonic isospin" - baryonic isospin density \n
 * \li "pion"             - pion density
 * \li "none"             - do not calculate density, print 0.0 \n
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

  const int ntest = config.take({"General", "Testparticles"}, 1);
  if (ntest <= 0) {
    throw std::invalid_argument("Invalid number of Testparticles "
                                "in config file!");
  }

  // The clock initializers are only read here and taken later when
  // assigning initial_clock_.
  return {{0.0f, config.read({"General", "Delta_Time"})},
          config.take({"Output", "Output_Interval"}), ntest,
          config.take({"General", "Gaussian_Sigma"}, 1.0),
          config.take({"General", "Gauss_Cutoff_In_Sigma"}, 4.0)};
}
}  // unnamed namespace

/**
 * Creates a verbose textual description of the setup of the Experiment.
 */
template <typename Modus>
std::ostream &operator<<(std::ostream &out, const Experiment<Modus> &e) {
  switch (e.time_step_mode_) {
    case TimeStepMode::None:
      out << "Not using time steps\n";
      break;
    case TimeStepMode::Fixed:
      out << "Using fixed time step size: "
          << e.parameters_.timestep_duration() << " fm/c\n";
      break;
  }
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
 * \key Force_Decays_At_End (bool, optional, default = true): \n
 * true - force all resonances to decay after last timestep \n
 * false - don't force decays (final output can contain resonances)
 *
 * \subpage pauliblocker
 */
template <typename Modus>
Experiment<Modus>::Experiment(Configuration config, bf::path output_path)
    : parameters_(create_experiment_parameters(config)),
      modus_(config["Modi"], parameters_),
      particles_(),
      nevents_(config.take({"General", "Nevents"})),
      end_time_(config.take({"General", "End_Time"})),
      delta_time_startup_(config.take({"General", "Delta_Time"})),
      force_decays_(
          config.take({"Collision_Term", "Force_Decays_At_End"}, true)),
      use_grid_(config.take({"General", "Use_Grid"}, true)),
      time_step_mode_(
          config.take({"General", "Time_Step_Mode"}, TimeStepMode::Fixed)) {
  const auto &log = logger<LogArea::Experiment>();
  log.info() << *this;

  // dilepton switch
  const bool dileptons_switch = config.take(
                                      {"Output", "Dileptons", "Enable"}, false);

  // create finders
  if (config.take({"Collision_Term", "Decays"}, true)) {
    action_finders_.emplace_back(new DecayActionsFinder());
  }
  if (dileptons_switch) {
    dilepton_finder_ = make_unique<DecayActionsFinderDilepton>();
  }
  if (config.take({"Collision_Term", "Collisions"}, true)) {
    action_finders_.emplace_back(new ScatterActionsFinder(config, parameters_));
  }
  if (config.has_value({"Collision_Term", "Pauli_Blocking"})) {
    log.info() << "Pauli blocking is ON.";
    pauli_blocker_ = make_unique<PauliBlocker>(
        config["Collision_Term"]["Pauli_Blocking"], parameters_);
  }

  // create outputs
  log.trace(source_location, " create OutputInterface objects");
  /*!\Userguide
    * \page input_output_options_ Output
    *
    * \key Output: \n
    * Below this key the configuration for the different output formats is
    * defined. To enable a certain output, set the 'Enable' key below the
    * selected format identifier. The identifiers are described below.
    * The following outputs exist:
    * \li \subpage input_oscar_particlelist
    * \li \subpage input_oscar_collisions
    * \li \subpage input_vtk
    * \li \subpage input_binary_collisions
    * \li \subpage input_binary_particles
    * \li \subpage input_root
    * \li \subpage input_dileptons
    */
  auto output_conf = config["Output"];
  /*!\Userguide
    * \page output_general_ Output files
    * There are different optional formats for SMASH output that are explained
    * below in more detail. Per default, the selected output files will be
    * saved in the directory ./data/\<run_id\>, where \<run_id\> is an integer
    * number starting from 0. At the beginning
    * of a run SMASH checks, if the ./data/0 directory exists. If it does not exist, it
    * is created and all output files are written there. If the directory
    * already exists, SMASH tries for ./data/1, ./data/2 and so on until it
    * finds a free number. The user can change output directory by a command
    * line option, if desired:
    * \code smash -o <user_output_dir> \endcode
    * SMASH supports several kinds of configurable output formats.
    * They are called OSCAR1999, OSCAR2013, binary OSCAR2013, VTK and ROOT
    * outputs. Every format can be switched on/off using option Enable in the
    * configuration file config.yaml. For more information on configuring the
    * output see corresponding pages: \ref input_oscar_particlelist,
    * \ref input_oscar_collisions, \ref input_binary_collisions,
    * \ref input_binary_particles, \ref input_root, \ref input_vtk.
    *
    * \key Details of output formats are explained here: \n
    * \li General block structure of OSCAR formats: \n
    *     \subpage oscar_general_
    * \li A family of OSCAR ASCII outputs.\n
    *     \subpage format_oscar_particlelist\n
    *     \subpage format_oscar_collisions
    * \li Binary outputs analoguous to OSCAR format\n
    *     \subpage format_binary_\n
    * \li Output in vtk format suitable for an easy
    *     visualization using paraview software:\n \subpage format_vtk
    * \li Formatted binary output that uses ROOT software
    *     (http://root.cern.ch).\n Fast to read and write, requires less
    *     disk space.\n \subpage format_root
    * \li Dilepton output in Oscar format: \n
    *     \subpage format_dilepton_output
    * \li \subpage collisions_output_in_box_modus_
    * \li \subpage output_vtk_lattice_
    */

  // loop until all OSCAR outputs are created (create_oscar_output will return
  // nullptr then).
  while (std::unique_ptr<OutputInterface> oscar =
              create_oscar_output(output_path, output_conf)) {
    outputs_.emplace_back(std::move(oscar));
  }
  if (static_cast<bool>(output_conf.take({"Vtk", "Enable"}))) {
    outputs_.emplace_back(new VtkOutput(output_path,
                                        std::move(output_conf["Vtk"])));
  } else {
    output_conf.take({"Vtk"});
  }
  if (static_cast<bool>(output_conf.take({"Binary_Collisions", "Enable"}))) {
    outputs_.emplace_back(new BinaryOutputCollisions(output_path,
                                  std::move(output_conf["Binary_Collisions"])));
  } else {
    output_conf.take({"Binary_Collisions"});
  }
  if (static_cast<bool>(output_conf.take({"Binary_Particles", "Enable"}))) {
    outputs_.emplace_back(new BinaryOutputParticles(output_path,
                                  std::move(output_conf["Binary_Particles"])));
  } else {
    output_conf.take({"Binary_Particles"});
  }
  if (static_cast<bool>(output_conf.take({"Root", "Enable"}))) {
#ifdef SMASH_USE_ROOT
    outputs_.emplace_back(new RootOutput(
                              output_path, std::move(output_conf["Root"])));
#else
    log.error() << "You requested Root output, but Root support has not been "
                    "compiled in.";
    output_conf.take({"Root"});
#endif
  } else {
    output_conf.take({"Root"});
  }
  if (static_cast<bool>(output_conf.take({"Density", "Enable"}))) {
    outputs_.emplace_back(new DensityOutput(output_path,
                              std::move(output_conf["Density"])));
  } else {
    output_conf.take({"Density"});
  }

  if (dileptons_switch) {
    dilepton_output_ = create_dilepton_output(output_path);
  }

  if (config.has_value({"Potentials"})) {
    if (time_step_mode_ == TimeStepMode::None) {
      log.error() << "Potentials only work with time steps!";
      throw std::invalid_argument("Can't use potentials without time steps!");
    }
    log.info() << "Potentials are ON.";
    // potentials need testparticles and gaussian sigma from parameters_
    potentials_ = make_unique<Potentials>(config["Potentials"], parameters_);
  }

  dens_type_ = config.take({"Output", "Density_Type"}, DensityType::None);
  log.info() << "Density type written to headers: " << dens_type_;

  /*!\Userguide
   * \page input_lattice_ Lattice
   *
   * \key Sizes (array<float,3>, required): \n
   *      Sizes of lattice in x, y, z directions in fm.
   *
   * \key Cell_Number (array<int,3>, required): \n
   *      Number of cells in x, y, z directions.
   *
   * \key Origin (array<float,3>, required): \n
   *      Coordinates of the left, down, near corner of the lattice in fm.
   *
   * \key Periodic (bool, required): \n
   *      Use periodic continuation or not. With periodic continuation
   *      x + i * lx is equivalent to x, same for y, z.
   *
   * \subpage input_vtk_lattice_
   *
   * For format of lattice output see \ref output_vtk_lattice_.
   *
   * \page input_vtk_lattice_ Printout
   *
   * User can print thermodynamical quantities on the lattice to vtk output.
   * For this one has to use the "Lattice: Printout" section of configuration.
   * Currently printing of custom density to vtk file is available.
   *
   * \key Density (string, optional, default = "none"): \n
   * Chooses which density to print.
   */

  // Create lattices
  if (config.has_value({"Lattice"})) {
    // Take lattice properties from config to assign them to all lattices
    const std::array<float, 3> l = config.take({"Lattice", "Sizes"});
    const std::array<int, 3> n = config.take({"Lattice", "Cell_Number"});
    const std::array<float, 3> origin = config.take({"Lattice", "Origin"});
    const bool periodic = config.take({"Lattice", "Periodic"});
    dens_type_lattice_printout_ = config.take(
                  {"Lattice", "Printout", "Density"}, DensityType::None);
    /* Create baryon and isospin density lattices regardless of config
       if potentials are on. This is because they allow to compute
       potentials faster */
    if (potentials_) {
      jmu_B_lat_ = make_unique<DensityLattice>(l, n, origin, periodic,
                                            LatticeUpdate::EveryTimestep);
      jmu_I3_lat_ = make_unique<DensityLattice>(l, n, origin, periodic,
                                              LatticeUpdate::EveryTimestep);
    } else {
      if (dens_type_lattice_printout_ == DensityType::Baryon) {
        jmu_B_lat_ = make_unique<DensityLattice>(l, n, origin, periodic,
                                                  LatticeUpdate::AtOutput);
      }
      if (dens_type_lattice_printout_ == DensityType::BaryonicIsospin) {
        jmu_I3_lat_ = make_unique<DensityLattice>(l, n, origin, periodic,
                                             LatticeUpdate::AtOutput);
      }
    }
    if (dens_type_lattice_printout_ != DensityType::None &&
        dens_type_lattice_printout_ != DensityType::BaryonicIsospin &&
        dens_type_lattice_printout_ != DensityType::Baryon) {
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
     << field<12, 3> << ((scatterings_total && time > really_small)
                          ? scatterings_total * 2 / (particles.size() * time)
                          : 0.)
     << field<10, 3> << scatterings_this_interval
     << field<12, 3> << particles.size()
     << field<10, 3> << elapsed_seconds;
  return ss.str();
}

template <typename Modus>
template <typename Container>
void Experiment<Modus>::perform_action(
    const ActionPtr &action, size_t &interactions_total,
    size_t &total_pauli_blocked, const Container &particles_before_actions) {
  const auto &log = logger<LogArea::Experiment>();
  if (action->is_valid(particles_)) {
    const ParticleList incoming_particles = action->incoming_particles();
    action->generate_final_state();
    ProcessType process_type = action->get_type();
    log.debug("Process Type is: ", process_type);
    if (pauli_blocker_ &&
        action->is_pauli_blocked(particles_, *pauli_blocker_.get())) {
      total_pauli_blocked++;
      return;
    }
    action->perform(&particles_, interactions_total);
    const ParticleList outgoing_particles = action->outgoing_particles();
    // Calculate Eckart rest frame density at the interaction point
    double rho = 0.0;
    if (dens_type_ != DensityType::None) {
      const FourVector r_interaction = action->get_interaction_point();
      constexpr bool compute_grad = false;
      rho = rho_eckart(r_interaction.threevec(), particles_before_actions,
                       parameters_, dens_type_, compute_grad).first;
    }
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

template <typename Modus>
void Experiment<Modus>::write_dilepton_action(const ActionPtr &action,
                                 const ParticleList &particles_before_actions) {
  if (action->is_valid(particles_)) {
    action->generate_final_state();
    // Calculate Eckart rest frame density at the interaction point
    const FourVector r_interaction = action->get_interaction_point();
    constexpr bool compute_grad = false;
    const double rho =
        rho_eckart(r_interaction.threevec(), particles_before_actions,
                   parameters_, dens_type_, compute_grad).first;
    // write dilepton output
    dilepton_output_->at_interaction(action->incoming_particles(),
                                     action->outgoing_particles(),
                                     rho,
                                     action->raw_weight_value(),
                                     action->get_type());
  }
}

template <typename Modus>
size_t Experiment<Modus>::run_time_evolution_without_time_steps(
    const int evt_num) {
  const auto &log = logger<LogArea::Experiment>();
  modus_.impose_boundary_conditions(&particles_);
  size_t interactions_total = 0, previous_interactions_total = 0,
         total_pauli_blocked = 0;

  // if no output is scheduled, trigger it manually
  if (!parameters_.is_output_time()) {
    log.info() << format_measurements(
        particles_, interactions_total, 0u,
        conserved_initial_, time_start_, parameters_.labclock.current_time());
  }

  const float start_time = parameters_.labclock.current_time();
  float time_left = end_time_ - start_time;

  // find actions for the initial list
  ParticleList search_list = particles_.copy_to_vector();
  Actions actions;
  for (const auto &finder : action_finders_) {
    actions.insert(finder->find_actions_in_cell(search_list, time_left));
  }

  // iterate over all actions
  while (!actions.is_empty()) {
    // get next action
    ActionPtr act = actions.pop();
    if (!act->is_valid(particles_)) {
      log.debug(~einhard::DRed(), "✘ ", act, " (discarded: invalid)");
      continue;
    }
    log.debug(~einhard::Green(), "✔ ", act);

    /* (1) Propagate to the next action. */

    const float action_time = act->time_of_execution();
    const float dt = action_time - parameters_.labclock.current_time();

    // we allow a very small negative time step that can result from imprecise
    // addition
    if (dt < -really_small) {
      log.error() << "dt = " << dt;
      throw std::runtime_error("Negative time step!");
    }

    float current_time;

    // only propagate the particles if dt is significantly larger than 0
    if (dt > really_small) {
      // set the time step according to our plan
      parameters_.labclock.set_timestep_duration(dt);

      // check if we need to do the intermediate output in the time until the
      // next action
      if (parameters_.need_intermediate_output()) {
        if (!parameters_.is_output_time()) {
          // we now set the clock to the output time and propagate the particles
          // until that time; then we do the output
          parameters_.set_timestep_for_next_output();
          ++parameters_.labclock;
          propagate_all();
          // after the output, the particles need to be propagated until the
          // action time
          const float remaining_dt =
              action_time - parameters_.labclock.current_time();
          parameters_.labclock.set_timestep_duration(remaining_dt);
        }
        intermediate_output(evt_num, interactions_total,
                            previous_interactions_total);
      }

      // set the clock manually instead of advancing it with the time step
      // to avoid loss of precision
      parameters_.labclock.reset(action_time);
      current_time = action_time;

      propagate_all();
    } else {
      // otherwise just keep the current time
      current_time = parameters_.labclock.current_time();
      parameters_.labclock.set_timestep_duration(0.f);
    }

    /* (2) Perform action. */

    // Update the positions of the incoming particles, because the information
    // in the action object will be outdated as the particles have been
    // propagated since the construction of the action.
    act->update_incoming(particles_);

    perform_action(act, interactions_total, total_pauli_blocked, particles_);
    modus_.impose_boundary_conditions(&particles_);

    /* (3) Check conservation laws. */

    std::string err_msg = conserved_initial_.report_deviations(particles_);
    if (!err_msg.empty()) {
      log.error() << err_msg;
      throw std::runtime_error("Violation of conserved quantities!");
    }

    /* (4) Find new actions. */

    time_left = end_time_ - current_time;
    const ParticleList& outgoing_particles = act->outgoing_particles();
    for (const auto &finder : action_finders_) {
      actions.insert(
          finder->find_actions_in_cell(outgoing_particles, time_left));
      actions.insert(finder->find_actions_with_surrounding_particles(
          outgoing_particles, particles_, time_left));
    }
  }
  // check if a final intermediate output is needed
  parameters_.labclock.end_tick_on_multiple(end_time_);
  ++parameters_.labclock;
  if (parameters_.is_output_time()) {
    propagate_all();
    intermediate_output(evt_num, interactions_total,
                        previous_interactions_total);
  }
  return interactions_total;
}

/* This is the loop over timesteps, carrying out collisions and decays
 * and propagating particles. */
template <typename Modus>
size_t Experiment<Modus>::run_time_evolution_fixed_time_step(
    const int evt_num) {
  const auto &log = logger<LogArea::Experiment>();
  modus_.impose_boundary_conditions(&particles_);
  size_t interactions_total = 0, previous_interactions_total = 0,
         total_pauli_blocked = 0;
  log.info() << format_measurements(
      particles_, interactions_total, 0u,
      conserved_initial_, time_start_, parameters_.labclock.current_time());

  Actions actions;
  Actions dilepton_actions;

  // minimal cell length for the grid
  const float min_cell_length = ScatterActionsFinder::min_cell_length(
      parameters_.testparticles, parameters_.timestep_duration());

  while (!(++parameters_.labclock > end_time_)) {
    /* (1.a) Create grid. */
    const auto &grid =
        use_grid_ ? modus_.create_grid(particles_, min_cell_length)
                  : modus_.create_grid(particles_, min_cell_length,
                                       CellSizeStrategy::Largest);
    /* (1.b) Iterate over cells and find actions. */
    grid.iterate_cells([&](const ParticleList &search_list) {
                         for (const auto &finder : action_finders_) {
                           actions.insert(finder->find_actions_in_cell(
                               search_list, parameters_.timestep_duration()));
                         }
                       },
                       [&](const ParticleList &search_list,
                           const ParticleList &neighbors_list) {
                         for (const auto &finder : action_finders_) {
                           actions.insert(finder->find_actions_with_neighbors(
                               search_list, neighbors_list,
                               parameters_.timestep_duration()));
                         }
                       });

    const auto particles_before_actions = particles_.copy_to_vector();

    /* (1.d) Dileptons */
    if (dilepton_finder_ != nullptr) {
      dilepton_actions.insert(dilepton_finder_->find_actions_in_cell(
                                              particles_before_actions,
                                              parameters_.timestep_duration()));

      if (!dilepton_actions.is_empty()) {
        while (!dilepton_actions.is_empty()) {
          write_dilepton_action(dilepton_actions.pop(),
                                particles_before_actions);
        }
      }
    }

    /* (2) Perform actions. */
    if (!actions.is_empty()) {
      while (!actions.is_empty()) {
        perform_action(actions.pop(), interactions_total, total_pauli_blocked,
                       particles_before_actions);
      }
      log.debug(~einhard::Blue(), particles_);
    } else {
      log.debug("no actions performed");
    }
    modus_.impose_boundary_conditions(&particles_);

    /* (3) Do propagation. */
    propagate_all();

    /* (4) Physics output during the run. */
    // if the timestep of labclock is different in the next tick than
    // in the current one, I assume it has been changed already. In that
    // case, I know what the next tick is and I can check whether the
    // output time is crossed within the next tick.
    if (parameters_.need_intermediate_output()) {
      intermediate_output(evt_num, interactions_total,
                          previous_interactions_total);
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
  return interactions_total;
}

template<typename Modus>
void Experiment<Modus>::intermediate_output(const int evt_num,
    size_t& interactions_total, size_t& previous_interactions_total) {
  const auto &log = logger<LogArea::Experiment>();
  const size_t interactions_this_interval =
      interactions_total - previous_interactions_total;
  previous_interactions_total = interactions_total;
  log.info() << format_measurements(
      particles_, interactions_total, interactions_this_interval,
      conserved_initial_, time_start_, parameters_.labclock.current_time());
  const LatticeUpdate lat_upd = LatticeUpdate::AtOutput;
  /* save evolution data */
  for (const auto &output : outputs_) {
    output->at_intermediate_time(particles_, evt_num, parameters_.labclock);
    // Thermodynamic output at some point versus time
    output->thermodynamics_output(particles_, parameters_);
    // Thermodynamic output on the lattice versus time
    switch (dens_type_lattice_printout_) {
      case DensityType::Baryon:
        update_density_lattice(jmu_B_lat_.get(), lat_upd,
                               DensityType::Baryon, parameters_, particles_);
        output->thermodynamics_output(std::string("rhoB"), *jmu_B_lat_,
                                                                 evt_num);
        break;
      case DensityType::BaryonicIsospin:
        update_density_lattice(jmu_I3_lat_.get(), lat_upd,
                     DensityType::BaryonicIsospin, parameters_, particles_);
        output->thermodynamics_output(std::string("rhoI3"), *jmu_I3_lat_,
                                                                 evt_num);
        break;
      case DensityType::None:
        break;
      default:
        update_density_lattice(jmu_custom_lat_.get(), lat_upd,
                       dens_type_lattice_printout_, parameters_, particles_);
        output->thermodynamics_output(std::string("rho"), *jmu_custom_lat_,
                                                                  evt_num);
    }
  }
}

template <typename Modus>
void Experiment<Modus>::propagate_all() {
  if (potentials_) {
    update_density_lattice(jmu_B_lat_.get(), LatticeUpdate::EveryTimestep,
                     DensityType::Baryon, parameters_, particles_);
    update_density_lattice(jmu_I3_lat_.get(), LatticeUpdate::EveryTimestep,
                     DensityType::BaryonicIsospin, parameters_, particles_);
    propagate(&particles_, parameters_, *potentials_);
  } else {
    propagate_straight_line(&particles_, parameters_);
  }
  modus_.impose_boundary_conditions(&particles_, outputs_);
}

template <typename Modus>
void Experiment<Modus>::do_final_decays(size_t &interactions_total) {
  size_t total_pauli_blocked = 0;

  // at end of time evolution: force all resonances to decay
  size_t interactions_old;
  do {
    Actions actions;
    Actions dilepton_actions;

    interactions_old = interactions_total;
    const auto particles_before_actions = particles_.copy_to_vector();

    /* Dileptons */
    if (dilepton_finder_ != nullptr) {
      dilepton_actions.insert(dilepton_finder_->find_final_actions(particles_));

      if (!dilepton_actions.is_empty()) {
        while (!dilepton_actions.is_empty()) {
          write_dilepton_action(dilepton_actions.pop(),
                                particles_before_actions);
        }
      }
    }
    /* Find actions. */
    for (const auto &finder : action_finders_) {
      actions.insert(finder->find_final_actions(particles_));
    }
    /* Perform actions. */
    while (!actions.is_empty()) {
      perform_action(actions.pop(), interactions_total, total_pauli_blocked,
                     particles_before_actions);
    }
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

template <typename Modus>
void Experiment<Modus>::final_output(size_t interactions_total,
                                     const int evt_num) {
  const auto &log = logger<LogArea::Experiment>();
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

  for (const auto &output : outputs_) {
    output->at_eventend(particles_, evt_num);
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
    size_t interactions_total;
    switch (time_step_mode_) {
      case TimeStepMode::None:
        interactions_total = run_time_evolution_without_time_steps(j);
        break;
      case TimeStepMode::Fixed:
        interactions_total = run_time_evolution_fixed_time_step(j);
        break;
    }

    if (force_decays_) {
      do_final_decays(interactions_total);
    }

    /* Output at event end */
    final_output(interactions_total, j);
  }
}

}  // namespace Smash
