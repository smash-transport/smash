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

  if (modus_chooser.compare("Box") == 0) {
    return make_unique<Experiment<BoxModus>>(config, output_path);
  } else if (modus_chooser.compare("List") == 0) {
    return make_unique<Experiment<ListModus>>(config, output_path);
  } else if (modus_chooser.compare("Collider") == 0) {
    return make_unique<Experiment<ColliderModus>>(config, output_path);
  } else if (modus_chooser.compare("Sphere") == 0) {
    return make_unique<Experiment<SphereModus>>(config, output_path);
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
 * files. To disable a certain output, comment the corresponding section out:
 *
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
    throw std::invalid_argument("Testparticle number should be positive!");
  }

  // If this Delta_Time option is absent (this can be for timestepless mode)
  // just assign 1.0 fm/c, reasonable value will be set at event initialization
  const double dt = config.take({"General", "Delta_Time"}, 1.0f);
  const double output_dt = config.take({"Output", "Output_Interval"});
  const bool two_to_one = config.take({"Collision_Term", "Two_to_One"}, true);
  const bool two_to_two = config.take({"Collision_Term", "Two_to_Two"}, true);
  const bool strings_switch = config.take({"Collision_Term", "Strings"}, false);
  const bool photons_switch = config.has_value({"Output", "Photons"}) ?
                    config.take({"Output", "Photons", "Enable"}, true) :
                    false;
  /// Elastic collisions between the nucleons with the square root s
  //  below low_snn_cut are excluded.
  const double low_snn_cut = config.take({"Collision_Term",
                                          "Elastic_NN_Cutoff_Sqrts"}, 1.98);
  const auto proton = ParticleType::try_find(pdg::p);
  const auto pion = ParticleType::try_find(pdg::pi_z);
  if (proton && pion &&
      low_snn_cut > proton->mass() + proton->mass() + pion->mass()) {
    log.warn("The cut-off should be below the threshold energy",
             " of the process: NN to NNpi");
  }
  return {{0.0f, dt}, {0.0, output_dt},
          ntest,
          config.take({"General", "Gaussian_Sigma"}, 1.0f),
          config.take({"General", "Gauss_Cutoff_In_Sigma"}, 4.0f),
          two_to_one,
          two_to_two,
          strings_switch,
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
          << e.parameters_.labclock.timestep_duration()
          << " fm/c\n";
      break;
    case TimeStepMode::Adaptive:
      out << "Using adaptive time steps, starting with: "
          << e.parameters_.labclock.timestep_duration()
          << " fm/c\n";
      break;
  }
  out << "End time: " << e.end_time_ << " fm/c\n";
  out << e.modus_;
  return out;
}

template <typename Modus>
template <typename TOutput>
void Experiment<Modus>::create_output(const char * name,
                   const bf::path &output_path,
                   Configuration&& conf) {
  const bool exists = conf.has_value_including_empty({name});
  if (!exists) {
    return;
  }
  if (conf.has_value({name, "Enable"})) {
    const auto &log = logger<LogArea::Experiment>();
    log.warn("Enable option is deprecated."
             " To disable/enable output comment/uncomment"
             " it out in the config.yaml.");
  }
  const bool enable = conf.take({name, "Enable"}, true);
  if (!enable) {
    conf.take({name});
    return;
  }
  outputs_.emplace_back(make_unique<TOutput>(output_path, conf[name]));
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
      delta_time_startup_(parameters_.labclock.timestep_duration()),
      force_decays_(
          config.take({"Collision_Term", "Force_Decays_At_End"}, true)),
      use_grid_(config.take({"General", "Use_Grid"}, true)),
      dileptons_switch_(config.has_value({"Output", "Dileptons"}) ?
                    config.take({"Output", "Dileptons", "Enable"}, true) :
                    false),
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
  if (parameters_.two_to_one || parameters_.two_to_two) {
    auto scat_finder = make_unique<ScatterActionsFinder>(config, parameters_,
                       nucleon_has_interacted_,
                       modus_.total_N_number(), modus_.proj_N_number(),
                       n_fractional_photons_);
    max_transverse_distance_sqr_ = scat_finder->max_transverse_distance_sqr(
                                                  parameters_.testparticles);
    action_finders_.emplace_back(std::move(scat_finder));
  }
  const float modus_l = modus_.length();
  if (modus_l > 0.f) {
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
    adaptive_parameters_ = make_unique<AdaptiveParameters>(
      config["General"]["Adaptive_Time_Step"]);
    log.info() << *adaptive_parameters_;
  }

  // create outputs
  log.trace(source_location, " create OutputInterface objects");

  auto output_conf = config["Output"];
  /*!\Userguide
    * \page output_general_ Output formats
    * Several different output formats are available in SMASH. They are
    * explained below in more detail. Per default, the selected output files
    * will be saved in the directory ./data/\<run_id\>, where \<run_id\> is an
    * integer number starting from 0. At the beginning of a run SMASH checks,
    * if the ./data/0 directory exists. If it does not exist, it is created and
    * all output files are written there. If the directory already exists,
    * SMASH tries for ./data/1, ./data/2 and so on until it finds a free
    * number. The user can change output directory by a command line option, if
    * desired:
    * \code smash -o <user_output_dir> \endcode
    * SMASH supports several kinds of configurable output formats.
    * They are called OSCAR1999, OSCAR2013, binary OSCAR2013, VTK and ROOT
    * outputs. Every format can be switched on/off by commenting/uncommenting
    * the corresponding section in the configuration file config.yaml. For more
    * information on configuring the output see corresponding pages: \ref
    * input_oscar_particlelist,
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
  while (OutputPtr oscar = create_oscar_output(output_path, output_conf)) {
    outputs_.emplace_back(std::move(oscar));
  }
  create_output<VtkOutput>("Vtk", output_path, std::move(output_conf));
  create_output<BinaryOutputCollisions>("Binary_Collisions",
                                        output_path, std::move(output_conf));
  create_output<BinaryOutputParticles>("Binary_Particles",
                                        output_path, std::move(output_conf));
#ifdef SMASH_USE_ROOT
  create_output<RootOutput>("Root", output_path, std::move(output_conf));
#else
  const bool enable_root = output_conf.take({"Root", "Enable"}, true);
  if (enable_root && output_conf.has_value_including_empty({"Root"})) {
    log.error("Root output requested, but Root support not compiled in");
  }
#endif
  create_output<ThermodynamicOutput>("Thermodynamics",
                                     output_path, std::move(output_conf));

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
   * written in the ouput in every single timestep without ever performing
   * them.  The are weighted with a "shining weight" to compensate for the
   * over-production.
   * \li The shining weight can be found in the weight element of the ouput.
   * \li The shining method is implemented in the DecayActionsFinderDilepton,
   * which is enabled together with the dilepton output.
   *
   * \note If you want dilepton decays, you also have to modify decaymodes.txt.
   * Dilepton decays are commented out by default.
   *
   * \key Format (string, required):\n
   * "Oscar" - The dilepton output is written to the file \c DileptonOutput.oscar
   * in \ref format_oscar_collisions (OSCAR2013 format) .\n
   * "Binary" - The dilepton output is written to the file \c DileptonOutput.bin
   * in \ref format_binary_ .\n
   * "Root" - The dilepton output is written to the file \c DileptonOutput.root
   * in \ref format_root .\n
   **/
  if (dileptons_switch_) {
    // create dilepton output object
    std::string format = config.take({"Output", "Dileptons", "Format"});
    if (format == "Oscar") {
      dilepton_output_ = create_dilepton_output(output_path);
    } else if (format == "Binary") {
      dilepton_output_ =
          make_unique<BinaryOutputCollisions>(output_path, "DileptonOutput");
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

  if (parameters_.photons_switch) {
    // create photon output object
    std::string format = config.take({"Output", "Photons", "Format"});
    if (format == "Oscar") {
      photon_output_ = create_photon_output(output_path);
    } else if (format == "Binary") {
      photon_output_ =
          make_unique<BinaryOutputCollisions>(output_path, "PhotonOutput");
    } else if (format == "Root") {
#ifdef SMASH_USE_ROOT
      photon_output_ = make_unique<RootOutput>(output_path, "PhotonOutput");
#else
      log.error() << "You requested Root output, but Root support has not been "
                     "compiled in.";
      output_conf.take({"Root"});
#endif
    } else {
      throw std::runtime_error("Bad Photon output format: " + format);
    }
  }

  // We can take away the Fermi motion flag, because the collider modus is
  // already initialized. We only need it when potentials are enabled, but we
  // always have to take it, otherwise SMASH will complain about unused
  // options.  We have to provide a default value for modi other than Collider.
  const FermiMotion motion = config.take({"Modi", "Collider", "Fermi_Motion"},
                                         FermiMotion::Off);
  if (config.has_value({"Potentials"})) {
    if (time_step_mode_ == TimeStepMode::None) {
      log.error() << "Potentials only work with time steps!";
      throw std::invalid_argument("Can't use potentials without time steps!");
    }
    if (motion == FermiMotion::Frozen) {
      log.error() << "Potentials don't work with frozen Fermi momenta! "
                     "Use normal Fermi motion instead.";
      throw std::invalid_argument("Can't use potentials "
                                  "with frozen Fermi momenta!");
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
    const std::array<float, 3> l = config.take({"Lattice", "Sizes"});
    const std::array<int, 3> n = config.take({"Lattice", "Cell_Number"});
    const std::array<float, 3> origin = config.take({"Lattice", "Origin"});
    const bool periodic = config.take({"Lattice", "Periodic"});
    dens_type_lattice_printout_ =
        config.take({"Lattice", "Printout", "Type"}, DensityType::None);
    const std::set<ThermodynamicQuantity> td_to_print =
        config.take({"Lattice", "Printout", "Quantities"});
    printout_tmn_ = (td_to_print.count(ThermodynamicQuantity::Tmn) > 0);
    printout_tmn_landau_ =
        (td_to_print.count(ThermodynamicQuantity::TmnLandau) > 0);
    printout_v_landau_ =
        (td_to_print.count(ThermodynamicQuantity::LandauVelocity) > 0);
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
      if (max_dt > 0.f && max_dt < timestep) {
        timestep = max_dt;
      }
      break;
  }
  Clock clock_for_this_event(start_time, timestep);
  parameters_.labclock = std::move(clock_for_this_event);

  /* Reset the output clock */
  const double dt_output = parameters_.outputclock.timestep_duration();
  const double zeroth_output_time = std::floor(start_time/dt_output)*dt_output;
  Clock output_clock(zeroth_output_time, dt_output);
  parameters_.outputclock = std::move(output_clock);

  log.debug("Lab clock: t_start = ", parameters_.labclock.current_time(),
           ", dt = ", parameters_.labclock.timestep_duration());
  log.debug("Output clock: t_start = ", parameters_.outputclock.current_time(),
           ", dt = ", parameters_.outputclock.timestep_duration());

  /* Save the initial conserved quantum numbers and total momentum in
   * the system for conservation checks */
  conserved_initial_ = QuantumNumbers(particles_);
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
                         ? 2.0 * scatterings_total / (particles.size() * time)
                         : 0.)
     << field<10, 3> << scatterings_this_interval
     << field<12, 3> << particles.size() << field<10, 3> << elapsed_seconds;
  return ss.str();
}

template <typename Modus>
template <typename Container>
bool Experiment<Modus>::perform_action(Action &action,
                                 const Container &particles_before_actions) {
  const auto &log = logger<LogArea::Experiment>();
  // Make sure to skip invalid and Pauli-blocked actions.
  if (!action.is_valid(particles_)) {
    log.debug(~einhard::DRed(), "✘ ", action, " (discarded: invalid)");
    return false;
  }
  action.generate_final_state();
  log.debug("Process Type is: ", action.get_type());
  if (pauli_blocker_ &&
      action.is_pauli_blocked(particles_, *pauli_blocker_)) {
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
    output->at_interaction(action, rho);
  }

  // At every collision photons can be produced.
  if (parameters_.photons_switch &&
      ScatterActionPhoton::is_photon_reaction(action.incoming_particles())) {
    // Time in the action constructor is relative to current time of incoming
    constexpr double action_time = 0.f;
    ScatterActionPhoton photon_act(action.incoming_particles(),
                                   action_time, n_fractional_photons_);
    // Add a completely dummy process to photon action.  The only important
    // thing is that its cross-section is equal to cross-section of action.
    // This can be done, because photon action is never performed, only
    // final state is generated and printed to photon output.
    photon_act.add_dummy_hadronic_channels(action.raw_weight_value());
    // Now add the actual photon reaction channel
    photon_act.add_single_channel();
    for (int i = 0; i < n_fractional_photons_; i++) {
      photon_act.generate_final_state();
      photon_output_->at_interaction(photon_act, rho);
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

  log.info() << format_measurements(particles_, interactions_total_, 0u,
                                    conserved_initial_, time_start_,
                                    parameters_.labclock.current_time());

  while (parameters_.labclock.current_time() < end_time_) {
    const double t = parameters_.labclock.current_time();
    const double dt = std::min(parameters_.labclock.timestep_duration(),
                              end_time_ - t);
    log.debug("Timestepless propagation for next ", dt, " fm/c.");

    /* (1.a) Create grid. */
    float min_cell_length = compute_min_cell_length(dt);
    log.debug("Creating grid with minimal cell length ", min_cell_length);
    const auto &grid = use_grid_
                           ? modus_.create_grid(particles_, min_cell_length)
                           : modus_.create_grid(particles_, min_cell_length,
                                                CellSizeStrategy::Largest);

    /* (1.b) Iterate over cells and find actions. */
    grid.iterate_cells(
        [&](const ParticleList &search_list) {
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

    /* (2) In case of adaptive timesteps adapt timestep size */
    if (time_step_mode_ ==  TimeStepMode::Adaptive && actions.size() > 0u) {
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

    ++parameters_.labclock;

    /* (5) Check conservation laws. */

    // Check conservation of conserved quantities if potentials and string
    // fragmentation are off.  If potentials are on then momentum is conserved
    // only in average.  If string fragmentation is on, then energy and
    // momentum are only very roughly conserved in high-energy collisions.
    if (!potentials_ && !parameters_.strings_switch) {
      std::string err_msg = conserved_initial_.report_deviations(particles_);
      if (!err_msg.empty()) {
        log.error() << err_msg;
        throw std::runtime_error("Violation of conserved quantities!");
      }
    }
  }

  if (pauli_blocker_) {
    log.info("Interactions: Pauli-blocked/performed = ", total_pauli_blocked_,
             "/", interactions_total_);
  }
}

template <typename Modus>
void Experiment<Modus>::propagate_and_shine(double to_time) {
  const double dt = propagate_straight_line(&particles_, to_time);
  if (dilepton_finder_ != nullptr) {
    dilepton_finder_->shine(particles_, dilepton_output_.get(), dt);
  }
}

template <typename Modus>
void Experiment<Modus>::run_time_evolution_timestepless(Actions& actions) {
  const auto &log = logger<LogArea::Experiment>();

  const double start_time = parameters_.labclock.current_time();
  const double end_time = std::min(parameters_.labclock.next_time(), end_time_);
  double time_left = end_time - start_time;
  log.debug("Timestepless propagation: ", "Actions size = ", actions.size(),
            ", start time = ", start_time,
            ", end time = ", end_time);

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
        log.debug(~einhard::DRed(), "✘ ", act, " (discarded: adaptive timestep"
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
    log.debug("Propagating until next action ", act, ", action time = ",
             act->time_of_execution());
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
  const uint64_t interactions_this_interval =
      interactions_total_ - previous_interactions_total_;
  previous_interactions_total_ = interactions_total_;
  log.info() << format_measurements(
      particles_, interactions_total_, interactions_this_interval,
      conserved_initial_, time_start_, parameters_.outputclock.current_time());
  const LatticeUpdate lat_upd = LatticeUpdate::AtOutput;
  /* save evolution data */
  for (const auto &output : outputs_) {
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
      dilepton_finder_->shine_final(particles_, dilepton_output_.get(), true);
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
    dilepton_finder_->shine_final(particles_, dilepton_output_.get(), false);
  }
}

template <typename Modus>
void Experiment<Modus>::final_output(const int evt_num) {
  const auto &log = logger<LogArea::Experiment>();
  // make sure the experiment actually ran (note: we should compare this
  // to the start time, but we don't know that. Therefore, we check that
  // the time is positive, which should heuristically be the same).
  if (likely(parameters_.labclock > 0)) {
    const uint64_t interactions_this_interval =
        interactions_total_ - previous_interactions_total_;
    log.info() << format_measurements(
      particles_, interactions_total_, interactions_this_interval,
      conserved_initial_, time_start_, parameters_.outputclock.current_time());
    log.info() << hline;
    log.info() << "Time real: " << SystemClock::now() - time_start_;
    /* if there are no particles no interactions happened */
    log.info() << "Final scattering rate: "
               << (particles_.is_empty() ? 0 : (2.0 * interactions_total_ /
                                                particles_.time() /
                                                particles_.size()))
               << " [fm-1]";
    log.info() << "Final interaction number: " << interactions_total_;
  }

  for (const auto &output : outputs_) {
    output->at_eventend(particles_, evt_num);
  }
  if (dilepton_output_ != nullptr) {
    dilepton_output_->at_eventend(particles_, evt_num);
  }
  if (photon_output_ != nullptr) {
    photon_output_->at_eventend(particles_, evt_num);
  }
}

template <typename Modus>
void Experiment<Modus>::run() {
  const auto &mainlog = logger<LogArea::Main>();
  for (int j = 0; j < nevents_; j++) {
    mainlog.info() << "Event " << j;

    /* Sample initial particles, start clock, some printout and book-keeping */
    initialize_new_event();
    /** In the ColliderMode, if the first collisions within the same nucleus are
     *  forbidden, then nucleon_has_interacted_ is created to record whether the nucleons inside
     *  the colliding nuclei have experienced any collisions or not */
    if (modus_.is_collider()) {
      if (!modus_.cll_in_nucleus()) {
        nucleon_has_interacted_.assign(modus_.total_N_number(), false);
      } else {
        nucleon_has_interacted_.assign(modus_.total_N_number(), true);
      }
    }
    /* Output at event start */
    for (const auto &output : outputs_) {
      output->at_eventstart(particles_, j);
    }
    if (dilepton_output_ != nullptr) {
      dilepton_output_->at_eventstart(particles_, j);
    }
    if (photon_output_ != nullptr) {
      photon_output_->at_eventstart(particles_, j);
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
