/*
 *
 *    Copyright (c) 2012-2017
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
#include "include/fourvector.h"
#include "include/listmodus.h"
#include "include/scatteractionphoton.h"
#include "include/scatteractionsfinder.h"
#include "include/spheremodus.h"
/* Outputs */
#include "include/binaryoutputcollisions.h"
#include "include/binaryoutputparticles.h"
#include "include/oscaroutput.h"
#include "include/thermodynamicoutput.h"
#ifdef SMASH_USE_ROOT
#include "include/rootoutput.h"
#endif
#include "include/vtkoutput.h"
#include "include/wallcrossingaction.h"

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
  using Seconds = chrono::duration<double>;
  using Minutes = chrono::duration<double, std::ratio<60>>;
  using Hours = chrono::duration<double, std::ratio<60 * 60>>;
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
ExperimentPtr ExperimentBase::create(Configuration config,
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

  if (modus_chooser == "Box") {
    return make_unique<Experiment<BoxModus>>(config, output_path);
  } else if (modus_chooser == "List") {
    return make_unique<Experiment<ListModus>>(config, output_path);
  } else if (modus_chooser == "Collider") {
    return make_unique<Experiment<ColliderModus>>(config, output_path);
  } else if (modus_chooser == "Sphere") {
    return make_unique<Experiment<SphereModus>>(config, output_path);
  } else {
    throw InvalidModusRequest("Invalid Modus (" + modus_chooser +
                              ") requested from ExperimentBase::create.");
  }
}

namespace {
/*!\Userguide
 * \page input_general_ General
 * \key Delta_Time (double, required): \n
 * Time step for the calculation, in fm/c.
 * Not required for timestepless mode.
 *
 * \key Testparticles (int, optional, default = 1): \n
 * How many test particles per real particles should be simulated.
 *
 * \key Gaussian_Sigma (double, optional, default 1.0): \n
 * Width [fm] of gaussians that represent Wigner density of particles.
 *
 * \key Gauss_Cutoff_In_Sigma (double, optional, default 4.0)
 * Distance in sigma at which gaussian is considered 0.
 *
 * \page input_output_options_ Output
 *
 * Description of options
 * ---------------------
 * \key Output_Interval (double, optional, default = End_Time): \n
 * Defines the period of intermediate output of the status of the simulated
 * system in Standard Output and other output formats which support this
 * functionality.
 *
 * \key Density_Type (string, optional, default = "none"): \n
 * Determines which kind of density is printed into the collision files.
 * Possible values:\n
 * \li "hadron"           - total hadronic density
 * \li "baryon"           - net baryon density
 * \li "baryonic isospin" - baryonic isospin density
 * \li "pion"             - pion density
 * \li "none"             - do not calculate density, print 0.0
 *
 * Futher options are defined for every single output \b content
 * (see \ref output_contents_ "output contents" for the list of
 * possible contents) in the following way:
 * \code
 * Content:
 *     Format: ["format1", "format2", ...]
 *     Option1: Value  # Content-specific options
 *     Option2: Value
 *     ...
 * \endcode
 *
 * To disable a certain output content,  remove or comment out the
 * corresponding section. Every output can be printed in several formats
 * simultaneously. The following option chooses list of formats:
 *
 * \key Format (list of formats, optional, default = []):\n
 * List of formats for writing particular content.
 * Possible formats for every content are listed and described in
 * \ref output_contents_ "output contents". List of available formats is
 * \ref list_of_output_formats "here".
 *
 * ### Content-specific output options
 * \anchor output_content_specific_options_
 *
 * - \b Particles
 *
 *   \key Extended (bool, optional, default = false): \n
 *   true - print extended information for each particle \n
 *   false - regular output for each particle
 *
 *   \key Only_Final (bool, optional, default = true): \n
 *   true - print only final particle list \n
 *   false - particle list at output interval including initial time
 * - \b Collisions
 *
 *   \key Extended (bool, optional, default = false): \n
 *   true - print extended information for each particle \n
 *   false - regular output for each particle
 *
 *   \key Print_Start_End (bool, optional, default = false): \n
 *   true - initial and final particle list is printed out \n
 *   false - initial and final particle list is not printed out
 * - \b Photons - see \ref input_photons
 * - \b Thermodynamics - see \subpage input_vtk_lattice_ for full spatial
 *   lattice output and \subpage ascii_thermodynamic_output_ for output at one
 *   point versus time.
 *
 * \anchor configuring_output_
 * Example configuring SMASH output
 * --------------
 * As an example, if one wants to have all of the following simultaneously:
 * \li particles at the end of event printed out in binary and Root formats
 * \li dileptons printed in Oscar2013 format
 * \li net baryon density at point (0, 0, 0) printed as a table
 *     against time every 1 fm/c
 *
 * then the output section of configuration will be the following.
 *
 * \code
 * Output:
 *     Output_Interval:  1.0
 *     Particles:
 *         Format:          ["Binary", "Root"]
 *         Only_Final:      True
 *     Dileptons:
 *         Format:          ["Oscar2013"]
 *     Thermodynamics:
 *         Format:          ["ASCII"]
 *         Type:            "baryon"
 *         Quantities:      ["rho_eckart"]
 *         Position:        [0.0, 0.0, 0.0]
 *         Smearing:        True
 * \endcode
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
    throw std::invalid_argument("Testparticle number should be positive!");
  }

  // If this Delta_Time option is absent (this can be for timestepless mode)
  // just assign 1.0 fm/c, reasonable value will be set at event initialization
  const double dt = config.take({"General", "Delta_Time"}, 1.);
  const double t_end = config.read({"General", "End_Time"});
  const double output_dt = config.take({"Output", "Output_Interval"}, t_end);
  const bool two_to_one = config.take({"Collision_Term", "Two_to_One"}, true);
  const bool two_to_two = config.take({"Collision_Term", "Two_to_Two"}, true);
  const bool strings_switch = config.take({"Collision_Term", "Strings"}, false);
  const NNbarTreatment nnbar_treatment = config.take(
      {"Collision_Term", "NNbar_Treatment"}, NNbarTreatment::NoAnnihilation);
  const bool photons_switch = config.has_value({"Output", "Photons"});
  /// Elastic collisions between the nucleons with the square root s
  //  below low_snn_cut are excluded.
  const double low_snn_cut =
      config.take({"Collision_Term", "Elastic_NN_Cutoff_Sqrts"}, 1.98);
  const auto proton = ParticleType::try_find(pdg::p);
  const auto pion = ParticleType::try_find(pdg::pi_z);
  if (proton && pion &&
      low_snn_cut > proton->mass() + proton->mass() + pion->mass()) {
    log.warn("The cut-off should be below the threshold energy",
             " of the process: NN to NNpi");
  }
  return {{0., dt},
          {0.0, output_dt},
          ntest,
          config.take({"General", "Gaussian_Sigma"}, 1.),
          config.take({"General", "Gauss_Cutoff_In_Sigma"}, 4.),
          two_to_one,
          two_to_two,
          strings_switch,
          nnbar_treatment,
          photons_switch,
          low_snn_cut};
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
          << e.parameters_.labclock.timestep_duration() << " fm/c\n";
      break;
    case TimeStepMode::Adaptive:
      out << "Using adaptive time steps, starting with: "
          << e.parameters_.labclock.timestep_duration() << " fm/c\n";
      break;
  }
  out << "End time: " << e.end_time_ << " fm/c\n";
  out << e.modus_;
  return out;
}

template <typename Modus>
void Experiment<Modus>::create_output(std::string format, std::string content,
                                      const bf::path &output_path,
                                      const OutputParameters &out_par) {
  const auto &log = logger<LogArea::Experiment>();
  log.info() << "Adding output " << content << " of format " << format
             << std::endl;

  if (format == "VTK" && content == "Particles") {
    outputs_.emplace_back(make_unique<VtkOutput>(output_path, content));
  } else if (format == "Root") {
#ifdef SMASH_USE_ROOT
    outputs_.emplace_back(
        make_unique<RootOutput>(output_path, content, out_par));
#else
    log.error("Root output requested, but Root support not compiled in");
#endif
  } else if (format == "Binary") {
    if (content == "Collisions" ||
        content == "Dileptons" || content == "Photons") {
      outputs_.emplace_back(
          make_unique<BinaryOutputCollisions>(output_path, content, out_par));
    } else if (content == "Particles") {
      outputs_.emplace_back(
          make_unique<BinaryOutputParticles>(output_path, content, out_par));
    }
  } else if (format == "Oscar1999" || format == "Oscar2013") {
    outputs_.emplace_back(
        create_oscar_output(format, content, output_path, out_par));
  } else if (content == "Thermodynamics" && format == "ASCII") {
    outputs_.emplace_back(
        make_unique<ThermodynamicOutput>(output_path, content, out_par));
  } else if (content == "Thermodynamics" && format == "VTK") {
    printout_lattice_td_ = true;
  } else {
    log.error() << "Unknown combination of format (" << format
                << ") and content (" << content << "). Fix the config.";
  }
}

/*!\Userguide
 * \page input_general_
 * \key End_Time (double, required): \n
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
 * \key Metric_Type (ExpansionMode, optional, default = NoExpansion): \n
 * NoExpansion - default SMASH run, with Minkowski metric \n
 * MasslessFRW - FRW expansion going as t^(1/2)
 * MassiveFRW - FRW expansion going as t^(2/3)
 * Exponential - FRW expansion going as e^(t/2)
 *
 * \key Expansion_Rate (double, optional, default = 0.1) \n
 * Corresponds to the speed of expansion of the universe in non minkowski
 * metrics \n
 * This value is useless if NoExpansion is selected; it corresponds to \n
 * \f$b_r/l_0\f$ if the metric type is MasslessFRW or MassiveFRW, and to \n
 * the parameter b in the Exponential expansion where \f$a(t) ~ e^{bt/2}\f$
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
      delta_time_startup_(parameters_.labclock.timestep_duration()),
      force_decays_(
          config.take({"Collision_Term", "Force_Decays_At_End"}, true)),
      use_grid_(config.take({"General", "Use_Grid"}, true)),
      metric_(
          config.take({"General", "Metric_Type"}, ExpansionMode::NoExpansion),
          config.take({"General", "Expansion_Rate"}, 0.1)),
      dileptons_switch_(config.has_value({"Output", "Dileptons"})),
      photons_switch_(config.has_value({"Output", "Photons"})),
      time_step_mode_(
          config.take({"General", "Time_Step_Mode"}, TimeStepMode::Fixed)) {
  const auto &log = logger<LogArea::Experiment>();
  log.info() << *this;

  // create finders
  if (dileptons_switch_) {
    dilepton_finder_ = make_unique<DecayActionsFinderDilepton>();
  }
  if (parameters_.photons_switch) {
    n_fractional_photons_ = config.take({"Output", "Photons", "Fractions"});
  }
  if (parameters_.two_to_one) {
    action_finders_.emplace_back(make_unique<DecayActionsFinder>());
  }
  if (parameters_.two_to_one || parameters_.two_to_two ||
      parameters_.strings_switch) {
    auto scat_finder = make_unique<ScatterActionsFinder>(
        config, parameters_, nucleon_has_interacted_, modus_.total_N_number(),
        modus_.proj_N_number(), n_fractional_photons_);
    max_transverse_distance_sqr_ =
        scat_finder->max_transverse_distance_sqr(parameters_.testparticles);
    action_finders_.emplace_back(std::move(scat_finder));
  } else {
    max_transverse_distance_sqr_ = maximum_cross_section / M_PI * fm2_mb;
  }
  const double modus_l = modus_.length();
  if (modus_l > 0.) {
    action_finders_.emplace_back(make_unique<WallCrossActionsFinder>(modus_l));
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
   * \key Smoothing_Factor (double, optional, default = 0.1) \n
   * Parameter of the exponential smoothing of the rate estimate.
   *
   * \key Target_Missed_Actions (double, optional, default = 0.01) \n
   * The fraction of missed actions that is targeted by the algorithm.
   *
   * \key Allowed_Deviation (double, optional, default = 2.5) \n
   * Limit by how much the target can be exceeded before the time step is
   * aborted.
   *
   **/
  if (time_step_mode_ == TimeStepMode::Adaptive) {
    adaptive_parameters_ = make_unique<AdaptiveParameters>(
        config["General"]["Adaptive_Time_Step"], delta_time_startup_);
    log.info() << *adaptive_parameters_;
  }

  // create outputs
  log.trace(source_location, " create OutputInterface objects");

  auto output_conf = config["Output"];
  /*!\Userguide
    * \page output_general_ Output
    *
    * Output directory
    * ----------------
    *
    * Per default, the selected output files
    * will be saved in the directory ./data/\<run_id\>, where \<run_id\> is an
    * integer number starting from 0. At the beginning of a run SMASH checks,
    * if the ./data/0 directory exists. If it does not exist, it is created and
    * all output files are written there. If the directory already exists,
    * SMASH tries for ./data/1, ./data/2 and so on until it finds a free
    * number.
    *
    * The user can change output directory by a command line option, if
    * desired:
    * \code smash -o <user_output_dir> \endcode
    *
    * Output content
    * --------------
    * \anchor output_contents_
    * Output in SMASH is distinguished by _content_ and _format_, where content
    * means the physical information contained in the output (e.g. list of
    * particles, list of interactions, thermodynamics, etc) and format (e.g.
    * Oscar, binary or ROOT). The same content can be printed out in several
    * formats _simultaneously_.
    *
    * For an example of choosing specific output contents see
    * \ref configuring_output_ "Configuring SMASH output".
    *
    * The list of possible contents follows:
    *
    * - \b Particles  List of particles at regular time intervals in the
    *                   computational frame or (optionally) only at the event end.
    *   - Available formats: \ref format_oscar_particlelist,
    *      \ref format_binary_, \ref format_root, \ref format_vtk
    * - \b Collisions List of interactions: collisions, decays, box wall
    *                 crossings and forced thermalizations. Information about
    *                 incoming, outgoing particles and the interaction itself
    *                 is printed out.
    *   - Available formats: \ref format_oscar_collisions, \ref format_binary_,
    *                 \ref format_root
    * - \b Dileptons  Special dilepton output, see \subpage input_dileptons.
    *   - Available formats: \ref format_oscar_collisions,
    *                   \ref format_binary_ and \ref format_root
    * - \b Photons    Special photon output, see \subpage input_photons.
    *   - Available formats: \ref format_oscar_collisions,
    *                   \ref format_binary_ and \ref format_root.
    * - \b Thermodynamics This output allows to print out thermodynamic
    *          quantities such as density, energy-momentum tensor,
    *          Landau velocity, etc at one selected point versus time
    *          (simple ASCII format table \subpage ascii_thermodynamic_output_)
    *          and on a spatial lattice  versus time (\ref output_vtk_lattice_).
    * \anchor list_of_output_formats
    * Output formats
    * --------------
    *
    * For choosing output formats see
    * \ref configuring_output_ "Configuring SMASH output".
    * Every output content can be printed out in several formats:
    * - \b "Oscar1999", \b "Oscar2013" - human-readable text output\n
    *   - For "Particles" content: \subpage format_oscar_particlelist
    *   - For "Collisions" content: \subpage format_oscar_collisions
    *   - General block structure of OSCAR formats: \subpage oscar_general_
    * - \b "Binary" - binary, not human-readable output
    *   - Faster to read and write than text outputs
    *   - Saves coordinates and momenta with the full double precision
    *   - General file structure is similar to \subpage oscar_general_
    *   - Detailed description: \subpage format_binary_
    * - \b "Root" - binary output in the format used by ROOT software
    *     (http://root.cern.ch)
    *   - Even faster to read and write, requires less disk space
    *   - Format description: \subpage format_root
    * - \b "VTK" - text output suitable for an easy
    *     visualization using paraview software
    *   - This output can be opened by paraview to see the visulalization.
    *   - For "Particles" content \subpage format_vtk
    *   - For "Thermodynamics" content \subpage output_vtk_lattice_
    * - \b "ASCII" - a human-readable text-format table of values
    *   - Used only for "Thermodynamics", see
    *     \subpage ascii_thermodynamic_output_
    *
    * \note Output of coordinates for the "Collisions" content in
    *       the periodic box has a feature:
    *       \subpage collisions_output_in_box_modus_
    */

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
   * written in the output at every hadron propagation without ever performing
   * them. The are weighted with a "shining weight" to compensate for the
   * over-production.
   * \li The shining weight can be found in the weight element of the output.
   * \li The shining method is implemented in the DecayActionsFinderDilepton,
   * which is enabled together with the dilepton output.
   *
   * \note If you want dilepton decays, you also have to modify decaymodes.txt.
   * Dilepton decays are commented out by default.
   **/

  /*!\Userguide
   * \page input_photons Photons
   * Todo(schaefer): document photons
   **/

  dens_type_ = config.take({"Output", "Density_Type"}, DensityType::None);
  log.info() << "Density type printed to headers: " << dens_type_;

  const OutputParameters output_parameters(std::move(output_conf));

  std::vector<std::string> output_contents = output_conf.list_upmost_nodes();
  for (const auto &content : output_contents) {
    auto this_output_conf = output_conf[content.c_str()];
    std::vector<std::string> formats = this_output_conf.take({"Format"});
    for (const auto &format : formats) {
      create_output(format, content, output_path, output_parameters);
    }
  }

  // We can take away the Fermi motion flag, because the collider modus is
  // already initialized. We only need it when potentials are enabled, but we
  // always have to take it, otherwise SMASH will complain about unused
  // options.  We have to provide a default value for modi other than Collider.
  const FermiMotion motion =
      config.take({"Modi", "Collider", "Fermi_Motion"}, FermiMotion::Off);
  if (config.has_value({"Potentials"})) {
    if (time_step_mode_ == TimeStepMode::None) {
      log.error() << "Potentials only work with time steps!";
      throw std::invalid_argument("Can't use potentials without time steps!");
    }
    if (motion == FermiMotion::Frozen) {
      log.error() << "Potentials don't work with frozen Fermi momenta! "
                     "Use normal Fermi motion instead.";
      throw std::invalid_argument(
          "Can't use potentials "
          "with frozen Fermi momenta!");
    }
    log.info() << "Potentials are ON.";
    // potentials need testparticles and gaussian sigma from parameters_
    potentials_ = make_unique<Potentials>(config["Potentials"], parameters_);
  }

  /*!\Userguide
   * \page input_lattice_ Lattice
   *
   * \key Sizes (array<double,3>, required): \n
   *      Sizes of lattice in x, y, z directions in fm.
   *
   * \key Cell_Number (array<int,3>, required): \n
   *      Number of cells in x, y, z directions.
   *
   * \key Origin (array<double,3>, required): \n
   *      Coordinates of the left, down, near corner of the lattice in fm.
   *
   * \key Periodic (bool, required): \n
   *      Use periodic continuation or not. With periodic continuation
   *      x + i * lx is equivalent to x, same for y, z.
   *
   * For format of lattice output see \ref output_vtk_lattice_. To configure
   * output of the quantities on the lattice to vtk files see
   * \ref input_output_options_.
   *
   * \page input_vtk_lattice_ lattice vtk output
   *
   * User can print thermodynamical quantities on the spatial lattice
   * to vtk output.
   * The lattice for the output is regulated by options of lattice
   * \subpage input_lattice_. The type of thermodynamic quantities is
   * chosen by the following options of the "Thermodynamic" output.
   *
   * \key Type (string, optional, default = "none"): \n
   * Chooses hadron/baryon/pion/baryonic isospin thermodynamic quantities
   *
   * \key Quantities (list of strings, optional, default = []): \n
   * List of quantities that can be printed:
   *  \li "rho_eckart": Eckart rest frame density
   *  \li "tmn": Energy-momentum tensor \f$T^{\mu\nu}(t,x,y,z) \f$
   *  \li "tmn_landau": Energy-momentum tensor in the Landau rest frame.
   *      This tensor is computed by boosting \f$T^{\mu\nu}(t,x,y,z) \f$
   *      to the local rest frame, where \f$T^{0i} \f$ = 0.
   *  \li "landau_velocity": Velocity of the Landau rest frame.
   *      The velocity is obtained from the energy-momentum tensor
   *      \f$T^{\mu\nu}(t,x,y,z) \f$ by solving the generalized eigenvalue
   *      equation \f$(T^{\mu\nu} - \lambda g^{\mu\nu})u_{\mu}=0 \f$.
  */

  // Create lattices
  if (config.has_value({"Lattice"})) {
    // Take lattice properties from config to assign them to all lattices
    const std::array<double, 3> l = config.take({"Lattice", "Sizes"});
    const std::array<int, 3> n = config.take({"Lattice", "Cell_Number"});
    const std::array<double, 3> origin = config.take({"Lattice", "Origin"});
    const bool periodic = config.take({"Lattice", "Periodic"});

    if (printout_lattice_td_) {
      dens_type_lattice_printout_ = output_parameters.td_dens_type;
      printout_tmn_ = output_parameters.td_tmn;
      printout_tmn_landau_ = output_parameters.td_tmn_landau;
      printout_v_landau_ = output_parameters.td_v_landau;
    }
    if (printout_tmn_ || printout_tmn_landau_ || printout_v_landau_) {
      Tmn_ = make_unique<RectangularLattice<EnergyMomentumTensor>>(
          l, n, origin, periodic, LatticeUpdate::AtOutput);
    }
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
      jmu_custom_lat_ = make_unique<DensityLattice>(l, n, origin, periodic,
                                                    LatticeUpdate::AtOutput);
    }
  } else if (printout_lattice_td_) {
    log.error("If you want Thermodynamic VTK output"
              ", configure a lattice for it.");
  }

  // Create forced thermalizer
  if (config.has_value({"Forced_Thermalization"})) {
    Configuration &&th_conf = config["Forced_Thermalization"];
    thermalizer_ = modus_.create_grandcan_thermalizer(th_conf);
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
  double start_time = modus_.initial_conditions(&particles_, parameters_);
  // For box modus make sure that particles are in the box. In principle, after
  // a correct initialization they should be, so this is just playing it safe.
  modus_.impose_boundary_conditions(&particles_, outputs_);
  /* Reset the simulation clock */
  double timestep = delta_time_startup_;

  switch (time_step_mode_) {
    case TimeStepMode::Fixed:
      break;
    case TimeStepMode::Adaptive:
      adaptive_parameters_->initialize(timestep);
      break;
    case TimeStepMode::None:
      timestep = end_time_ - start_time;
      // Take care of the box modus + timestepless propagation
      const double max_dt = modus_.max_timestep(max_transverse_distance_sqr_);
      if (max_dt > 0. && max_dt < timestep) {
        timestep = max_dt;
      }
      break;
  }
  Clock clock_for_this_event(start_time, timestep);
  parameters_.labclock = std::move(clock_for_this_event);

  /* Reset the output clock */
  const double dt_output = parameters_.outputclock.timestep_duration();
  const double zeroth_output_time =
      std::floor(start_time / dt_output) * dt_output;
  Clock output_clock(zeroth_output_time, dt_output);
  parameters_.outputclock = std::move(output_clock);

  log.debug("Lab clock: t_start = ", parameters_.labclock.current_time(),
            ", dt = ", parameters_.labclock.timestep_duration());
  log.debug("Output clock: t_start = ", parameters_.outputclock.current_time(),
            ", dt = ", parameters_.outputclock.timestep_duration());

  /* Save the initial conserved quantum numbers and total momentum in
   * the system for conservation checks */
  conserved_initial_ = QuantumNumbers(particles_);
  wall_actions_total_ = 0;
  previous_wall_actions_total_ = 0;
  interactions_total_ = 0;
  previous_interactions_total_ = 0;
  total_pauli_blocked_ = 0;
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
  ss << field<5> << time << field<12, 3> << difference.momentum().x0()
     << field<12, 3> << difference.momentum().abs3()
     << field<12, 3> << (time > really_small
                             ? 2.0 * scatterings_total /
                                   (particles.size() * time)
                             : 0.)
     << field<10, 3> << scatterings_this_interval
     << field<12, 3> << particles.size() << field<10, 3> << elapsed_seconds;
  return ss.str();
}

template <typename Modus>
template <typename Container>
bool Experiment<Modus>::perform_action(
    Action &action, const Container &particles_before_actions) {
  const auto &log = logger<LogArea::Experiment>();
  // Make sure to skip invalid and Pauli-blocked actions.
  if (!action.is_valid(particles_)) {
    log.debug(~einhard::DRed(), "✘ ", action, " (discarded: invalid)");
    return false;
  }
  action.generate_final_state();
  log.debug("Process Type is: ", action.get_type());
  if (pauli_blocker_ && action.is_pauli_blocked(particles_, *pauli_blocker_)) {
    total_pauli_blocked_++;
    return false;
  }
  if (modus_.is_collider()) {
    // Mark incoming nucleons as interacted - now they are permitted
    // to collide with nucleons from their native nucleus
    for (const auto &incoming : action.incoming_particles()) {
      assert(incoming.id() >= 0);
      if (incoming.id() < modus_.total_N_number()) {
        nucleon_has_interacted_[incoming.id()] = true;
      }
    }
  }
  // Make sure to pick a non-zero integer, because 0 is reserved for "no
  // interaction yet".
  const auto id_process = static_cast<uint32_t>(interactions_total_ + 1);
  action.perform(&particles_, id_process);
  interactions_total_++;
  if (action.get_type() == ProcessType::Wall) {
    wall_actions_total_++;
  }
  // Calculate Eckart rest frame density at the interaction point
  double rho = 0.0;
  if (dens_type_ != DensityType::None) {
    const FourVector r_interaction = action.get_interaction_point();
    constexpr bool compute_grad = false;
    rho = rho_eckart(r_interaction.threevec(), particles_before_actions,
                     density_param_, dens_type_, compute_grad)
              .first;
  }
  /*!\Userguide
   * \page collisions_output_in_box_modus_ Collision output in box modus
   * \note When SMASH is running in the box modus, particle coordinates
   * in the collision output can be out of the box. This is not an error.  Box
   * boundary conditions are intentionally not imposed before collision output
   * to allow unambiguous finding of the interaction point.
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
    if (!output->is_dilepton_output() && !output->is_photon_output()) {
      output->at_interaction(action, rho);
    }
  }

  // At every collision photons can be produced.
  if (photons_switch_ &&
      (ScatterActionPhoton::is_photon_reaction(action.incoming_particles()) !=
       ScatterActionPhoton::ReactionType::no_reaction)) {
    // Time in the action constructor is relative to
    // current time of incoming
    constexpr double action_time = 0.;
    ScatterActionPhoton photon_act(action.incoming_particles(), action_time,
                                   n_fractional_photons_);
    // Add a completely dummy process to photon action.  The only important
    // thing is that its cross-section is equal to cross-section of action.
    // This can be done, because photon action is never performed, only
    // final state is generated and printed to photon output.
    photon_act.add_dummy_hadronic_channels(action.raw_weight_value());
    // Now add the actual photon reaction channel
    photon_act.add_single_channel();
    for (int i = 0; i < n_fractional_photons_; i++) {
      photon_act.generate_final_state();
      for (const auto &output : outputs_) {
        if (output->is_photon_output()) {
          output->at_interaction(photon_act, rho);
        }
      }
    }
  }
  log.debug(~einhard::Green(), "✔ ", action);
  return true;
}

/// Make sure `interactions_total` can be represented as a 32-bit integer.
/// This is necessary for converting to a `id_process`. The latter is 32-bit
/// integer, because it is written like this to binary output.
static void check_interactions_total(uint64_t interactions_total) {
  constexpr uint64_t max_uint32 = std::numeric_limits<uint32_t>::max();
  if (interactions_total >= max_uint32) {
    throw std::runtime_error("Integer overflow in total interaction number!");
  }
}

template <typename Modus>
void Experiment<Modus>::run_time_evolution() {
  Actions actions;

  const auto &log = logger<LogArea::Experiment>();
  const auto &log_ad_ts = logger<LogArea::AdaptiveTS>();

  log.info() << format_measurements(
      particles_, interactions_total_ - wall_actions_total_, 0u,
      conserved_initial_, time_start_, parameters_.labclock.current_time());

  while (parameters_.labclock.current_time() < end_time_) {
    const double t = parameters_.labclock.current_time();
    const double dt =
        std::min(parameters_.labclock.timestep_duration(), end_time_ - t);
    log.debug("Timestepless propagation for next ", dt, " fm/c.");

    /* Perform forced thermalization if required */
    if (thermalizer_ &&
        thermalizer_->is_time_to_thermalize(parameters_.labclock)) {
      const bool ignore_cells_under_treshold = true;
      thermalizer_->update_lattice(particles_, density_param_,
                                   ignore_cells_under_treshold);
      const double current_t = parameters_.labclock.current_time();
      thermalizer_->thermalize(particles_, current_t,
                               parameters_.testparticles);
      ThermalizationAction th_act(*thermalizer_, current_t);
      if (th_act.any_particles_thermalized()) {
        perform_action(th_act, particles_);
      }
    }

    /* (1.a) Create grid. */
    double min_cell_length = compute_min_cell_length(dt);
    log.debug("Creating grid with minimal cell length ", min_cell_length);
    const auto &grid = use_grid_
                           ? modus_.create_grid(particles_, min_cell_length)
                           : modus_.create_grid(particles_, min_cell_length,
                                                CellSizeStrategy::Largest);

    /* (1.b) Iterate over cells and find actions. */
    grid.iterate_cells(
        [&](const ParticleList &search_list) {
          for (const auto &finder : action_finders_) {
            actions.insert(finder->find_actions_in_cell(search_list, dt));
          }
        },
        [&](const ParticleList &search_list,
            const ParticleList &neighbors_list) {
          for (const auto &finder : action_finders_) {
            actions.insert(finder->find_actions_with_neighbors(
                search_list, neighbors_list, dt));
          }
        });

    /* (2) In case of adaptive timesteps adapt timestep size */
    if (time_step_mode_ == TimeStepMode::Adaptive && actions.size() > 0u) {
      double new_timestep = parameters_.labclock.timestep_duration();
      if (adaptive_parameters_->update_timestep(actions, particles_.size(),
                                                &new_timestep)) {
        parameters_.labclock.set_timestep_duration(new_timestep);
        log_ad_ts.info("New timestep is set to ", new_timestep);
      }
    }

    /* (3) Propagation from action to action until the end of timestep */
    run_time_evolution_timestepless(actions);

    /* (4) Update potentials (if computed on the lattice) and
           compute new momenta according to equations of motion */
    if (potentials_) {
      update_potentials();
      update_momenta(&particles_, parameters_.labclock.timestep_duration(),
                     *potentials_, dUB_dr_lat_.get(), dUI3_dr_lat_.get());
    }

    /* (5) Expand universe if non-minkowskian metric; updates
           positions and momenta according to the selected expansion */
    if (metric_.mode_ != ExpansionMode::NoExpansion) {
      expand_space_time(&particles_, parameters_, metric_);
    }

    ++parameters_.labclock;

    /* (5) Check conservation laws. */

    // Check conservation of conserved quantities if potentials and string
    // fragmentation are off.  If potentials are on then momentum is conserved
    // only in average.  If string fragmentation is on, then energy and
    // momentum are only very roughly conserved in high-energy collisions.
    if (!potentials_ && !parameters_.strings_switch &&
        metric_.mode_ == ExpansionMode::NoExpansion) {
      std::string err_msg = conserved_initial_.report_deviations(particles_);
      if (!err_msg.empty()) {
        log.error() << err_msg;
        throw std::runtime_error("Violation of conserved quantities!");
      }
    }
  }

  if (pauli_blocker_) {
    log.info("Interactions: Pauli-blocked/performed = ", total_pauli_blocked_,
             "/", interactions_total_ - wall_actions_total_);
  }
}

template <typename Modus>
void Experiment<Modus>::propagate_and_shine(double to_time) {
  const double dt =
      propagate_straight_line(&particles_, to_time, beam_momentum_);
  if (dilepton_finder_ != nullptr) {
    for (const auto &output : outputs_) {
      dilepton_finder_->shine(particles_, output.get(), dt);
    }
  }
}

template <typename Modus>
void Experiment<Modus>::run_time_evolution_timestepless(Actions &actions) {
  const auto &log = logger<LogArea::Experiment>();

  const double start_time = parameters_.labclock.current_time();
  const double end_time = std::min(parameters_.labclock.next_time(), end_time_);
  double time_left = end_time - start_time;
  log.debug("Timestepless propagation: ", "Actions size = ", actions.size(),
            ", start time = ", start_time, ", end time = ", end_time);

  // iterate over all actions
  while (!actions.is_empty()) {
    // get next action
    ActionPtr act = actions.pop();
    if (!act->is_valid(particles_)) {
      log.debug(~einhard::DRed(), "✘ ", act, " (discarded: invalid)");
      continue;
    }
    if (act->time_of_execution() > end_time) {
      if (time_step_mode_ == TimeStepMode::Adaptive) {
        log.debug(~einhard::DRed(), "✘ ", act,
                  " (discarded: adaptive timestep"
                  " mode decreased timestep and this action is too late)");
      } else {
        log.error(act, " scheduled later than end time: t_action[fm/c] = ",
                  act->time_of_execution(), ", t_end[fm/c] = ", end_time);
      }
    }
    log.debug(~einhard::Green(), "✔ ", act);

    while (next_output_time() <= act->time_of_execution()) {
      log.debug("Propagating until output time: ", next_output_time());
      propagate_and_shine(next_output_time());
      ++parameters_.outputclock;
      intermediate_output();
    }

    /* (1) Propagate to the next action. */
    log.debug("Propagating until next action ", act,
              ", action time = ", act->time_of_execution());
    propagate_and_shine(act->time_of_execution());

    /* (2) Perform action. */

    // Update the positions of the incoming particles, because the information
    // in the action object will be outdated as the particles have been
    // propagated since the construction of the action.
    act->update_incoming(particles_);

    const bool performed = perform_action(*act, particles_);

    // No need to update actions for outgoing particles
    // if the action is not performed.
    if (!performed) {
      continue;
    }
    const auto particles_before_actions = particles_.copy_to_vector();

    /* (3) Update actions for newly-produced particles. */

    time_left = end_time - act->time_of_execution();
    const ParticleList &outgoing_particles = act->outgoing_particles();
    for (const auto &finder : action_finders_) {
      // Outgoing particles can still decay, cross walls...
      actions.insert(
          finder->find_actions_in_cell(outgoing_particles, time_left));
      // ... and collide with other particles.
      actions.insert(finder->find_actions_with_surrounding_particles(
          outgoing_particles, particles_, time_left));
    }

    check_interactions_total(interactions_total_);
  }

  while (next_output_time() <= end_time) {
    log.debug("Propagating until output time: ", next_output_time());
    propagate_and_shine(next_output_time());
    ++parameters_.outputclock;
    // Avoid duplicating printout at event end time
    if (parameters_.outputclock.current_time() < end_time_) {
      intermediate_output();
    }
  }

  log.debug("Propagating to time ", end_time);
  propagate_and_shine(end_time);
}

template <typename Modus>
void Experiment<Modus>::intermediate_output() {
  const auto &log = logger<LogArea::Experiment>();
  const uint64_t wall_actions_this_interval =
      wall_actions_total_ - previous_wall_actions_total_;
  previous_wall_actions_total_ = wall_actions_total_;
  const uint64_t interactions_this_interval = interactions_total_ -
                                              previous_interactions_total_ -
                                              wall_actions_this_interval;
  previous_interactions_total_ = interactions_total_;
  log.info() << format_measurements(
      particles_, interactions_total_ - wall_actions_total_,
      interactions_this_interval, conserved_initial_, time_start_,
      parameters_.outputclock.current_time());
  const LatticeUpdate lat_upd = LatticeUpdate::AtOutput;
  /*if (thermalizer_) {
    thermalizer_->update_lattice(particles_, density_param_);
    thermalizer_->print_statistics(parameters_.labclock);
  }*/
  /* save evolution data */
  for (const auto &output : outputs_) {
    if (output->is_dilepton_output() || output->is_photon_output()) {
      continue;
    }
    output->at_intermediate_time(particles_, parameters_.outputclock,
                                 density_param_);

    // Thermodynamic output on the lattice versus time
    switch (dens_type_lattice_printout_) {
      case DensityType::Baryon:
        update_density_lattice(jmu_B_lat_.get(), lat_upd, DensityType::Baryon,
                               density_param_, particles_);
        output->thermodynamics_output(ThermodynamicQuantity::EckartDensity,
                                      DensityType::Baryon, *jmu_B_lat_);
        break;
      case DensityType::BaryonicIsospin:
        update_density_lattice(jmu_I3_lat_.get(), lat_upd,
                               DensityType::BaryonicIsospin, density_param_,
                               particles_);
        output->thermodynamics_output(ThermodynamicQuantity::EckartDensity,
                                      DensityType::BaryonicIsospin,
                                      *jmu_I3_lat_);
        break;
      case DensityType::None:
        break;
      default:
        update_density_lattice(jmu_custom_lat_.get(), lat_upd,
                               dens_type_lattice_printout_, density_param_,
                               particles_);
        output->thermodynamics_output(ThermodynamicQuantity::EckartDensity,
                                      dens_type_lattice_printout_,
                                      *jmu_custom_lat_);
    }
    if (printout_tmn_ || printout_tmn_landau_ || printout_v_landau_) {
      update_Tmn_lattice(Tmn_.get(), lat_upd, dens_type_lattice_printout_,
                         density_param_, particles_);
      if (printout_tmn_) {
        output->thermodynamics_output(ThermodynamicQuantity::Tmn,
                                      dens_type_lattice_printout_, *Tmn_);
      }
      if (printout_tmn_landau_) {
        output->thermodynamics_output(ThermodynamicQuantity::TmnLandau,
                                      dens_type_lattice_printout_, *Tmn_);
      }
      if (printout_v_landau_) {
        output->thermodynamics_output(ThermodynamicQuantity::LandauVelocity,
                                      dens_type_lattice_printout_, *Tmn_);
      }
    }

    if (thermalizer_) {
      output->thermodynamics_output(*thermalizer_);
    }
  }
}

template <typename Modus>
void Experiment<Modus>::update_potentials() {
  if (potentials_) {
    if (potentials_->use_skyrme() && jmu_B_lat_ != nullptr) {
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
                             DensityType::BaryonicIsospin, density_param_,
                             particles_);
      const size_t UI3lattice_size = UI3_lat_->size();
      for (size_t i = 0; i < UI3lattice_size; i++) {
        (*UI3_lat_)[i] = potentials_->symmetry_pot((*jmu_I3_lat_)[i].density());
      }
      UI3_lat_->compute_gradient_lattice(dUI3_dr_lat_.get());
    }
  }
}

template <typename Modus>
void Experiment<Modus>::do_final_decays() {
  /* At end of time evolution: Force all resonances to decay. In order to handle
   * decay chains, we need to loop until no further actions occur. */
  uint64_t interactions_old;
  const auto particles_before_actions = particles_.copy_to_vector();
  do {
    Actions actions;

    interactions_old = interactions_total_;

    /* Dileptons: shining of remaining resonances */
    if (dilepton_finder_ != nullptr) {
      for (const auto &output : outputs_) {
        dilepton_finder_->shine_final(particles_, output.get(), true);
      }
    }
    /* Find actions. */
    for (const auto &finder : action_finders_) {
      actions.insert(finder->find_final_actions(particles_));
    }
    /* Perform actions. */
    while (!actions.is_empty()) {
      perform_action(*actions.pop(), particles_before_actions);
    }
    // loop until no more decays occur
  } while (interactions_total_ > interactions_old);

  /* Dileptons: shining of stable particles at the end */
  if (dilepton_finder_ != nullptr) {
    for (const auto &output : outputs_) {
      dilepton_finder_->shine_final(particles_, output.get(), false);
    }
  }
}

template <typename Modus>
void Experiment<Modus>::final_output(const int evt_num) {
  const auto &log = logger<LogArea::Experiment>();
  // make sure the experiment actually ran (note: we should compare this
  // to the start time, but we don't know that. Therefore, we check that
  // the time is positive, which should heuristically be the same).
  if (likely(parameters_.labclock > 0)) {
    const uint64_t wall_actions_this_interval =
        wall_actions_total_ - previous_wall_actions_total_;
    const uint64_t interactions_this_interval = interactions_total_ -
                                                previous_interactions_total_ -
                                                wall_actions_this_interval;
    log.info() << format_measurements(
        particles_, interactions_total_ - wall_actions_total_,
        interactions_this_interval, conserved_initial_, time_start_,
        parameters_.outputclock.current_time());
    log.info() << hline;
    log.info() << "Time real: " << SystemClock::now() - time_start_;
    /* if there are no particles no interactions happened */
    log.info() << "Final scattering rate: "
               << (particles_.is_empty()
                       ? 0
                       : (2.0 * (interactions_total_ - wall_actions_total_) /
                          particles_.time() / particles_.size()))
               << " [fm-1]";
    log.info() << "Final interaction number: "
               << interactions_total_ - wall_actions_total_;
    // Check if there are unformed particles
    int unformed_particles_count = 0;
    for (const auto & particle : particles_) {
      if (particle.formation_time() > end_time_) {
        unformed_particles_count++;
      }
    }
    if (unformed_particles_count > 0) {
      log.warn("End time might be too small. ", unformed_particles_count,
               " unformed particles were found at the end of the evolution.");
    }
  }

  for (const auto &output : outputs_) {
    output->at_eventend(particles_, evt_num, modus_.impact_parameter());
  }
}

template <typename Modus>
void Experiment<Modus>::run() {
  const auto &mainlog = logger<LogArea::Main>();
  for (int j = 0; j < nevents_; j++) {
    mainlog.info() << "Event " << j;

    /* Sample initial particles, start clock, some printout and book-keeping */
    initialize_new_event();
    /* In the ColliderModus, if the first collisions within the same nucleus are
     * forbidden, then nucleon_has_interacted_ is created to record whether the
     * nucleons inside
     * the colliding nuclei have experienced any collisions or not */
    if (modus_.is_collider()) {
      if (!modus_.cll_in_nucleus()) {
        nucleon_has_interacted_.assign(modus_.total_N_number(), false);
      } else {
        nucleon_has_interacted_.assign(modus_.total_N_number(), true);
      }
    }
    /* In the ColliderModus, if Fermi motion is frozen, assign the beam momenta
     * to
     * the nucleons in both the projectile and the target. */
    if (modus_.is_collider() && modus_.fermi_motion() == FermiMotion::Frozen) {
      for (int i = 0; i < modus_.total_N_number(); i++) {
        const auto mass_beam = particles_.copy_to_vector()[i].effective_mass();
        const auto v_beam = i < modus_.proj_N_number()
                                ? modus_.velocity_projectile()
                                : modus_.velocity_target();
        const auto gamma = 1.0 / std::sqrt(1.0 - v_beam * v_beam);
        beam_momentum_.emplace_back(FourVector(gamma * mass_beam, 0.0, 0.0,
                                               gamma * v_beam * mass_beam));
      }
    }

    /* Output at event start */
    for (const auto &output : outputs_) {
      output->at_eventstart(particles_, j);
    }

    run_time_evolution();

    if (force_decays_) {
      do_final_decays();
    }

    /* Output at event end */
    final_output(j);
  }
}

}  // namespace Smash
