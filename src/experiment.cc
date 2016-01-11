/*
 *
 *    Copyright (c) 2012-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/experiment.h"

#include "include/actions.h"
#include "include/boxmodus.h"
#include "include/collidermodus.h"
#include "include/cxx14compat.h"
#include "include/decayactionsfinder.h"
#include "include/decayactionsfinderdilepton.h"
#include "include/listmodus.h"
#include "include/propagation.h"
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
                                                  const bf::path &output_path) {
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
 * Not required for timestepless mode.
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
 * \li "hadron"           - total hadronic density
 * \li "baryon"           - net baryon density
 * \li "baryonic isospin" - baryonic isospin density
 * \li "pion"             - pion density
 * \li "none"             - do not calculate density, print 0.0
 *
 * The output section has several subsections, relating to different output
 * files. To enable a certain output, set the 'Enable' key in the corresponding
 * subsection:
 * \li \subpage input_oscar_particlelist
 * \li \subpage input_oscar_collisions
 * \li \subpage input_vtk
 * \li \subpage input_binary_collisions
 * \li \subpage input_binary_particles
 * \li \subpage input_root
 * \li \subpage input_dileptons
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

  const float dt = (config.has_value({"General", "Time_Step_Mode"}) &&
             config.read({"General", "Time_Step_Mode"}) == TimeStepMode::None)
             ? 0.0f
             : config.take({"General", "Delta_Time"});
  return {{0.0f, dt},
          config.take({"Output", "Output_Interval"}), ntest,
          config.take({"General", "Gaussian_Sigma"}, 1.0f),
          config.take({"General", "Gauss_Cutoff_In_Sigma"}, 4.0f)};
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
    case TimeStepMode::Adaptive:
      out << "Using adaptive time steps, starting with: "
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
 * \key Use_Grid (bool, optional, default = true): \n
 * true - a grid is used to reduce the combinatorics of interaction lookup \n
 * false - no grid is used
 *
 * \key Time_Step_Mode (string, optional, default = Fixed): \n
 * The mode of time stepping. Possible values: \n
 * None - No time steps are used. Cannot be used with potentials \n
 * Fixed - Fixed-sized time steps \n
 * Adaptive - Time steps with adaptive sizes
 *
 * \page input_collision_term_ Collision_Term
 *
 * \key Two_to_One (bool, optional, default = true) \n
 * Enable 2 <--> 1 processes (resonance formation and decays).
 *
 * \key Two_to_Two (bool, optional, default = true) \n
 * Enable 2 <--> 2 collisions.
 *
 * \key Force_Decays_At_End (bool, optional, default = true): \n
 * true - force all resonances to decay after last timestep \n
 * false - don't force decays (final output can contain resonances)
 *
 * \subpage pauliblocker
 */
template <typename Modus>
Experiment<Modus>::Experiment(Configuration config, const bf::path &output_path)
    : parameters_(create_experiment_parameters(config)),
      density_param_(DensityParameters(parameters_)),
      modus_(config["Modi"], parameters_),
      particles_(),
      nevents_(config.take({"General", "Nevents"})),
      end_time_(config.take({"General", "End_Time"})),
      delta_time_startup_(parameters_.timestep_duration()),
      force_decays_(
          config.take({"Collision_Term", "Force_Decays_At_End"}, true)),
      use_grid_(config.take({"General", "Use_Grid"}, true)),
      time_step_mode_(
          config.take({"General", "Time_Step_Mode"}, TimeStepMode::Fixed)) {
  const auto &log = logger<LogArea::Experiment>();
  log.info() << *this;

  const bool two_to_one = config.take({"Collision_Term", "Two_to_One"}, true);
  const bool two_to_two = config.take({"Collision_Term", "Two_to_Two"}, true);
  const bool strings_switch = config.take({"Collision_Term", "Strings"}, true);
  const bool dileptons_switch = config.take(
                                      {"Output", "Dileptons", "Enable"}, false);

  // create finders
  if (two_to_one) {
    action_finders_.emplace_back(new DecayActionsFinder());
  }
  if (two_to_one || two_to_two) {
    auto scat_finder = make_unique<ScatterActionsFinder>(config, parameters_,
                                                       two_to_one, two_to_two,
                                                       strings_switch);
    max_transverse_distance_sqr_ = scat_finder->max_transverse_distance_sqr(
                                                    parameters_.testparticles);
    action_finders_.emplace_back(std::move(scat_finder));
  }
  if (dileptons_switch) {
    dilepton_finder_ = make_unique<DecayActionsFinderDilepton>();
  }
  if (config.has_value({"Collision_Term", "Pauli_Blocking"})) {
    log.info() << "Pauli blocking is ON.";
    pauli_blocker_ = make_unique<PauliBlocker>(
        config["Collision_Term"]["Pauli_Blocking"], parameters_);
  }

  /*!\Userguide
   * \page input_general_ General
   * \subpage input_general_adaptive_
   * (optional)
   *
   * \page input_general_adaptive_ Adaptive_Time_Step
   * Additional parameters for the adaptive time step mode.
   *
   * \key Smoothing_Factor (float, optional, default = 0.1) \n
   * Parameter of the exponential smoothing of the rate estimate.
   *
   * \key Target_Missed_Actions (float, optional, default = 0.01) \n
   * The fraction of missed actions that is targeted by the algorithm.
   *
   * \key Allowed_Deviation (float, optional, default = 2.5) \n
   * Limit by how much the target can be exceeded before the time step is
   * aborted.
   *
   **/
  if (time_step_mode_ == TimeStepMode::Adaptive) {
    std::unique_ptr<AdaptiveParameters> adapt_params =
        make_unique<AdaptiveParameters>();
    if (config.has_value(
            {"General", "Adaptive_Time_Step", "Smoothing_Factor"})) {
      adapt_params->smoothing_factor =
          config.take({"General", "Adaptive_Time_Step", "Smoothing_Factor"});
    }
    if (config.has_value(
            {"General", "Adaptive_Time_Step", "Target_Missed_Actions"})) {
      adapt_params->target_missed_actions = config.take(
          {"General", "Adaptive_Time_Step", "Target_Missed_Actions"});
    }
    if (config.has_value(
            {"General", "Adaptive_Time_Step", "Allowed_Deviation"})) {
      adapt_params->deviation_factor =
          config.take({"General", "Adaptive_Time_Step", "Allowed_Deviation"});
    }
    log.info("Parameters for the adaptive time step:\n",
             "  Smoothing factor: ", adapt_params->smoothing_factor, "\n",
             "  Target missed actions: ",
             100 * adapt_params->target_missed_actions, "%", "\n",
             "  Allowed deviation: ", adapt_params->deviation_factor);
    adaptive_parameters_ = std::move(adapt_params);
  }

  // create outputs
  log.trace(source_location, " create OutputInterface objects");

  auto output_conf = config["Output"];
  /*!\Userguide
    * \page output_general_ Output formats
    * Several different output formats are available in SMASH. They are explained
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

  /*!\Userguide
   * \page input_dileptons Dileptons
   * Enables Dilepton Output together with DecayActionsFinderDilepton.
   * Dilepton Output saves information about decays, which include Dileptons,
   * at every timestep.
   *
   * The treatment of Dilepton Decays is special:
   *
   * \li Dileptons are treated via the time integration method, also called
   * 'shining', as described in \iref{Schmidt:2008hm}, chapter 2D.
   * This means that, because dilepton decays are so rare, possible decays are
   * written in the ouput in every single timestep without ever performing them.
   * The are weighted with a "shining weight" to compensate for the over-production.
   * \li The shining weight can be found in the weight element of the ouput.
   * \li The shining method is implemented in the DecayActionsFinderDilepton,
   * which is enabled together with the dilepton output.
   *
   * \note If you want dilepton decays, you also have to modify decaymodes.txt.
   * Dilepton decays are commented out by default.
   *
   * \key Enable (bool, optional, default = false):\n
   * true - Dilepton Output and DecayActionsFinderDilepton enabled\n
   * false - no Dilepton Output and no DecayActionsFinderDilepton
   *
   * \key Format (string, required):\n
   * "Oscar" - The dilepton output is written to the file \c DileptonOutput.oscar
   * in \ref format_oscar_collisions (OSCAR2013 format) .\n
   * "Binary" - The dilepton output is written to the file \c DileptonOutput.bin
   * in \ref format_binary_ .\n
   * "Root" - The dilepton output is written to the file \c DileptonOutput.root
   * in \ref format_root .\n
   **/
  if (dileptons_switch) {
    // create dilepton output object
    std::string format = config.take({"Output", "Dileptons", "Format"});
    if (format == "Oscar") {
      dilepton_output_ = create_dilepton_output(output_path);
    } else if (format == "Binary") {
      dilepton_output_ = make_unique<BinaryOutputCollisions>(output_path,
                                                             "DileptonOutput");
    } else if (format == "Root") {
#ifdef SMASH_USE_ROOT
      dilepton_output_ = make_unique<RootOutput>(output_path, "DileptonOutput");
#else
    log.error() << "You requested Root output, but Root support has not been "
                   "compiled in.";
    output_conf.take({"Root"});
#endif
    } else {
        throw std::runtime_error("Bad dilepton output format: " + format);
    }
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
      if (potentials_->use_skyrme()) {
        jmu_B_lat_ = make_unique<DensityLattice>(l, n, origin, periodic,
                                              LatticeUpdate::EveryTimestep);
        UB_lat_ = make_unique<RectangularLattice<double>>(
                       l, n, origin, periodic, LatticeUpdate::EveryTimestep);
        dUB_dr_lat_ = make_unique<RectangularLattice<ThreeVector>>(
                       l, n, origin, periodic, LatticeUpdate::EveryTimestep);
      }
      if (potentials_->use_symmetry()) {
        jmu_I3_lat_ = make_unique<DensityLattice>(l, n, origin, periodic,
                                              LatticeUpdate::EveryTimestep);
        UI3_lat_ = make_unique<RectangularLattice<double>>(
                        l, n, origin, periodic, LatticeUpdate::EveryTimestep);
        dUI3_dr_lat_ = make_unique<RectangularLattice<ThreeVector>>(
                        l, n, origin, periodic, LatticeUpdate::EveryTimestep);
      }
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
                                       uint64_t scatterings_total,
                                       uint64_t scatterings_this_interval,
                                       const QuantumNumbers &conserved_initial,
                                       SystemTimePoint time_start,
                                       double time) {
  const SystemTimeSpan elapsed_seconds = SystemClock::now() - time_start;

  const QuantumNumbers current_values(particles);
  const QuantumNumbers difference = conserved_initial - current_values;

  std::ostringstream ss;
  ss << field<5> << time
     << field<12, 3> << difference.momentum().x0()
     << field<12, 3> << difference.momentum().abs3()
     << field<12, 3> << (time > really_small
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
    Action &action, uint64_t &interactions_total,
    uint64_t &total_pauli_blocked, const Container &particles_before_actions) {
  const auto &log = logger<LogArea::Experiment>();
  if (!action.is_valid(particles_)) {
    log.debug(~einhard::DRed(), "✘ ", action, " (discarded: invalid)");
    return;
  }
  action.generate_final_state();
  log.debug("Process Type is: ", action.get_type());
  if (pauli_blocker_ &&
      action.is_pauli_blocked(particles_, *pauli_blocker_.get())) {
    total_pauli_blocked++;
    return;
  }
  // Make sure to pick a non-zero integer, because 0 is reserved for "no
  // interaction yet".
  const auto id_process = static_cast<uint32_t>(interactions_total + 1);
  action.perform(&particles_, id_process);
  interactions_total++;
  // Calculate Eckart rest frame density at the interaction point
  double rho = 0.0;
  if (dens_type_ != DensityType::None) {
    const FourVector r_interaction = action.get_interaction_point();
    constexpr bool compute_grad = false;
    rho = rho_eckart(r_interaction.threevec(), particles_before_actions,
                     density_param_, dens_type_, compute_grad).first;
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
    output->at_interaction(action, rho);
  }
  log.debug(~einhard::Green(), "✔ ", action);
}

template <typename Modus>
void Experiment<Modus>::write_dilepton_action(Action &action,
                                 const ParticleList &particles_before_actions) {
  if (action.is_valid(particles_)) {
    action.generate_final_state();
    // Calculate Eckart rest frame density at the interaction point
    const FourVector r_interaction = action.get_interaction_point();
    constexpr bool compute_grad = false;
    const double rho =
        rho_eckart(r_interaction.threevec(), particles_before_actions,
                   density_param_, dens_type_, compute_grad).first;
    // write dilepton output
    dilepton_output_->at_interaction(action, rho);
  }
}

/// Make sure `interactions_total` can be represented as a 32-bit integer.
/// This is necessary for converting to a `id_process`.
static void check_interactions_total(uint64_t interactions_total) {
  constexpr uint64_t max_uint32 = std::numeric_limits<uint32_t>::max();
  if (interactions_total >= max_uint32) {
    throw std::runtime_error("Integer overflow in total interaction number!");
  }
}

template <typename Modus>
uint64_t Experiment<Modus>::run_time_evolution_without_time_steps() {
  const auto &log = logger<LogArea::Experiment>();
  modus_.impose_boundary_conditions(&particles_);
  uint64_t interactions_total = 0, previous_interactions_total = 0,
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
        intermediate_output(interactions_total, previous_interactions_total);
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

    perform_action(*act, interactions_total, total_pauli_blocked, particles_);
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

    check_interactions_total(interactions_total);
  }
  // check if a final intermediate output is needed
  parameters_.labclock.end_tick_on_multiple(end_time_);
  ++parameters_.labclock;
  if (parameters_.is_output_time()) {
    propagate_all();
    intermediate_output(interactions_total, previous_interactions_total);
  }
  return interactions_total;
}

/* This is the loop over timesteps, carrying out collisions and decays
 * and propagating particles. */
template <typename Modus>
uint64_t Experiment<Modus>::run_time_evolution_fixed_time_step() {
  const auto &log = logger<LogArea::Experiment>();
  modus_.impose_boundary_conditions(&particles_);
  uint64_t interactions_total = 0, previous_interactions_total = 0,
         total_pauli_blocked = 0;
  log.info() << format_measurements(
      particles_, interactions_total, 0u,
      conserved_initial_, time_start_, parameters_.labclock.current_time());

  Actions actions;
  Actions dilepton_actions;
  const float dt = parameters_.timestep_duration();
  // minimal cell length of the grid for collision finding
  const float min_cell_length = compute_min_cell_length(dt);

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
                               search_list, dt));
                         }
                       },
                       [&](const ParticleList &search_list,
                           const ParticleList &neighbors_list) {
                         for (const auto &finder : action_finders_) {
                           actions.insert(finder->find_actions_with_neighbors(
                               search_list, neighbors_list, dt));
                         }
                       });

    const auto particles_before_actions = particles_.copy_to_vector();

    /* (1.d) Dileptons */
    if (dilepton_finder_ != nullptr) {
      dilepton_actions.insert(dilepton_finder_->find_actions_in_cell(
                                              particles_before_actions, dt));

      if (!dilepton_actions.is_empty()) {
        while (!dilepton_actions.is_empty()) {
          write_dilepton_action(*dilepton_actions.pop(),
                                particles_before_actions);
        }
      }
    }

    /* (2) Perform actions. */
    if (!actions.is_empty()) {
      while (!actions.is_empty()) {
        perform_action(*actions.pop(), interactions_total, total_pauli_blocked,
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
    if (parameters_.need_intermediate_output()) {
      intermediate_output(interactions_total, previous_interactions_total);
    }
    // Check conservation of conserved quantities if potentials are off.
    // If potentials are on then momentum is conserved only in average

    #ifdef PYTHIA_FOUND
    if (!potentials_) {
      std::string err_msg = conserved_initial_.report_deviations(particles_);
      if (!err_msg.empty()) {
        log.error() << err_msg;
        throw std::runtime_error("Violation of conserved quantities!");
      }
    }
    #endif
    check_interactions_total(interactions_total);
  }

  if (pauli_blocker_) {
    log.info("Interactions: Pauli-blocked/performed = ", total_pauli_blocked,
             "/", interactions_total);
  }
  return interactions_total;
}

/* This is the loop over timesteps, carrying out collisions and decays
 * and propagating particles. */
template <typename Modus>
uint64_t Experiment<Modus>::run_time_evolution_adaptive_time_steps(
                                const AdaptiveParameters &adaptive_parameters) {
  const auto &log = logger<LogArea::Experiment>();
  const auto &log_ad_ts = logger<LogArea::AdaptiveTS>();
  modus_.impose_boundary_conditions(&particles_);
  uint64_t interactions_total = 0, previous_interactions_total = 0,
           total_pauli_blocked = 0;
  bool observed_first_action = false;

  // if there is no output scheduled at the beginning, trigger it manually
  if (!parameters_.is_output_time()) {
    log.info() << format_measurements(
        particles_, interactions_total, 0u,
        conserved_initial_, time_start_, parameters_.labclock.current_time());
  }

  float rate =
      adaptive_parameters.rate_from_dt(parameters_.timestep_duration());
  Actions actions;
  uint32_t num_time_steps = 0u;
  float min_dt = std::numeric_limits<float>::infinity();
  while (parameters_.labclock.current_time() < end_time_) {
    num_time_steps++;
    if (parameters_.labclock.next_time() > end_time_) {
      // set the time step size such that we stop at end_time_
      parameters_.labclock.end_tick_on_multiple(end_time_);
    }
    float dt = parameters_.timestep_duration();

    /* (1.a) Create grid. */
    float min_cell_length = compute_min_cell_length(dt);
    const auto &grid = modus_.create_grid(particles_, min_cell_length);

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

    /* (2) Calculate time step size. */
    log_ad_ts.debug() << hline;
    if (!observed_first_action && actions.size() > 0u) {
      log_ad_ts.debug("First interaction.");
      observed_first_action = true;
    }

    bool end_early = false;
    if (observed_first_action) {
      float fraction_missed;
      float allowed_deviation;
      std::tie(fraction_missed, allowed_deviation) =
          adaptive_parameters.calc_missed_actions_allowed_deviation(
              actions, rate, particles_.size());
      const float current_rate = fraction_missed / dt;
      const float rate_deviation = current_rate - rate;
      // check if the current rate deviates too strongly from the expected value
      if (rate_deviation > allowed_deviation) {
        log_ad_ts.debug("End time step early.");
        end_early = true;
        parameters_.labclock.set_timestep_duration(
            adaptive_parameters.new_dt(current_rate));
      }
      // update the estimate of the rate
      rate += adaptive_parameters.smoothing_factor * rate_deviation;
      // set the size of the next time step
      dt = adaptive_parameters.new_dt(rate);
    }

    if (!observed_first_action) {
      log_ad_ts.debug("Averaged rate: ", 0.f);
    } else {
      log_ad_ts.debug("Averaged rate: ", rate);
    }
    log_ad_ts.debug("Time step size: ",
                    parameters_.labclock.timestep_duration());
    const float this_dt = parameters_.labclock.timestep_duration();
    if (this_dt < min_dt) {
      min_dt = this_dt;
    }

    /* (3) Physics output during the run. */
    if (parameters_.need_intermediate_output()) {
      // The following block is necessary because we want to make sure that the
      // intermediate output happens at exactly the requested time and not
      // slightly before.
      if (!parameters_.is_output_time()) {
        const float next_time = parameters_.labclock.next_time();
        // set the time step such that it ends on the next output time
        parameters_.set_timestep_for_next_output();
        ++parameters_.labclock;

        // perform actions until the output time
        const auto particles_before_actions = particles_.copy_to_vector();
        while (!actions.is_empty()) {
          auto action = actions.pop();
          if (action->time_of_execution() >
              parameters_.labclock.current_time()) {
            // reinsert action
            actions.insert(std::move(action));
            break;
          }
          perform_action(*action, interactions_total, total_pauli_blocked,
                         particles_before_actions);
        }
        modus_.impose_boundary_conditions(&particles_);
        propagate_all();
        for (const ActionPtr &action : actions) {
          if (action->is_valid(particles_)) {
            action->update_incoming(particles_);
          }
        }
        parameters_.labclock.end_tick_on_multiple(next_time);
      }

      intermediate_output(interactions_total, previous_interactions_total);
    }

    ++parameters_.labclock;

    /* (4) Perform actions. */
    if (!actions.is_empty()) {
      const auto particles_before_actions = particles_.copy_to_vector();
      while (!actions.is_empty()) {
        auto action = actions.pop();
        if (end_early &&
            action->time_of_execution() > parameters_.labclock.current_time()) {
          actions.clear();
          log_ad_ts.debug("Actions discarded because of early ending.");
          break;
        }
        perform_action(*action, interactions_total, total_pauli_blocked,
                       particles_before_actions);
      }
      log.debug(~einhard::Blue(), particles_);
    } else {
      log.debug("no actions performed");
    }
    modus_.impose_boundary_conditions(&particles_);

    /* (5) Do propagation. */
    propagate_all();

    /* (6) Set duration of next time step. */
    parameters_.labclock.set_timestep_duration(dt);

    // Check conservation of conserved quantities if potentials are off.
    // If potentials are on then momentum is conserved only in average
    if (!potentials_) {
      std::string err_msg = conserved_initial_.report_deviations(particles_);
      if (!err_msg.empty()) {
        log.error() << err_msg;
        throw std::runtime_error("Violation of conserved quantities!");
      }
    }

    check_interactions_total(interactions_total);
  }

  // check if a final intermediate output is needed
  if (parameters_.is_output_time()) {
    intermediate_output(interactions_total, previous_interactions_total);
  }

  if (pauli_blocker_) {
    log.info("Collisions: pauliblocked/total = ", total_pauli_blocked, "/",
             interactions_total);
  }
  log.info("Number of time steps = ", num_time_steps);
  log.info("Smallest time step size = ", min_dt);
  return interactions_total;
}

template<typename Modus>
void Experiment<Modus>::intermediate_output(uint64_t& interactions_total,
                                        uint64_t& previous_interactions_total) {
  const auto &log = logger<LogArea::Experiment>();
  const uint64_t interactions_this_interval =
      interactions_total - previous_interactions_total;
  previous_interactions_total = interactions_total;
  log.info() << format_measurements(
      particles_, interactions_total, interactions_this_interval,
      conserved_initial_, time_start_, parameters_.labclock.current_time());
  const LatticeUpdate lat_upd = LatticeUpdate::AtOutput;
  /* save evolution data */
  for (const auto &output : outputs_) {
    output->at_intermediate_time(particles_, parameters_.labclock,
                                 density_param_);

    // Thermodynamic output on the lattice versus time
    switch (dens_type_lattice_printout_) {
      case DensityType::Baryon:
        update_density_lattice(jmu_B_lat_.get(), lat_upd,
                               DensityType::Baryon, density_param_, particles_);
        output->thermodynamics_output(std::string("rhoB"), *jmu_B_lat_);
        break;
      case DensityType::BaryonicIsospin:
        update_density_lattice(jmu_I3_lat_.get(), lat_upd,
                     DensityType::BaryonicIsospin, density_param_, particles_);
        output->thermodynamics_output(std::string("rhoI3"), *jmu_I3_lat_);
        break;
      case DensityType::None:
        break;
      default:
        update_density_lattice(jmu_custom_lat_.get(), lat_upd,
                       dens_type_lattice_printout_, density_param_, particles_);
        output->thermodynamics_output(std::string("rho"), *jmu_custom_lat_);
    }
  }
}

template <typename Modus>
void Experiment<Modus>::propagate_all() {
  if (potentials_) {
    if (potentials_->use_skyrme() && jmu_B_lat_!= nullptr) {
      update_density_lattice(jmu_B_lat_.get(), LatticeUpdate::EveryTimestep,
                       DensityType::Baryon, density_param_, particles_);
      const size_t UBlattice_size = UB_lat_->size();
      for (size_t i = 0; i < UBlattice_size; i++) {
        (*UB_lat_)[i] = potentials_->skyrme_pot((*jmu_B_lat_)[i].density());
      }
      UB_lat_->compute_gradient_lattice(dUB_dr_lat_.get());
    }
    if (potentials_->use_symmetry() && jmu_I3_lat_ != nullptr) {
      update_density_lattice(jmu_I3_lat_.get(), LatticeUpdate::EveryTimestep,
                      DensityType::BaryonicIsospin, density_param_, particles_);
      const size_t UI3lattice_size = UI3_lat_->size();
      for (size_t i = 0; i < UI3lattice_size; i++) {
        (*UI3_lat_)[i] = potentials_->symmetry_pot(
                                        (*jmu_I3_lat_)[i].density());
      }
      UI3_lat_->compute_gradient_lattice(dUI3_dr_lat_.get());
    }
    propagate(&particles_, parameters_, *potentials_,
              dUB_dr_lat_.get(), dUI3_dr_lat_.get());
  } else {
    propagate_straight_line(&particles_, parameters_);
  }
  modus_.impose_boundary_conditions(&particles_, outputs_);
}

template <typename Modus>
void Experiment<Modus>::do_final_decays(uint64_t &interactions_total) {
  uint64_t total_pauli_blocked = 0;

  // at end of time evolution: force all resonances to decay
  uint64_t interactions_old;
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
          write_dilepton_action(*dilepton_actions.pop(),
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
      perform_action(*actions.pop(), interactions_total, total_pauli_blocked,
                     particles_before_actions);
    }
    // loop until no more decays occur
  } while (interactions_total > interactions_old);

  /* Do one final propagation step. */
  if (potentials_) {
    propagate(&particles_, parameters_, *potentials_,
              dUB_dr_lat_.get(), dUI3_dr_lat_.get());
  } else {
    propagate_straight_line(&particles_, parameters_);
  }
  modus_.impose_boundary_conditions(&particles_, outputs_);
}

template <typename Modus>
void Experiment<Modus>::final_output(uint64_t interactions_total,
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
    log.info() << "Final interaction number: " << interactions_total;
  }

  for (const auto &output : outputs_) {
    output->at_eventend(particles_, evt_num);
  }
  if (dilepton_output_ != nullptr) {
    dilepton_output_->at_eventend(particles_, evt_num);
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
    if (dilepton_output_ != nullptr) {
      dilepton_output_->at_eventstart(particles_, j);
    }

    /* the time evolution of the relevant subsystem */
    uint64_t interactions_total;
    switch (time_step_mode_) {
      case TimeStepMode::None:
        interactions_total = run_time_evolution_without_time_steps();
        break;
      case TimeStepMode::Fixed:
        interactions_total = run_time_evolution_fixed_time_step();
        break;
      case TimeStepMode::Adaptive:
        interactions_total =
            run_time_evolution_adaptive_time_steps(*adaptive_parameters_);
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
