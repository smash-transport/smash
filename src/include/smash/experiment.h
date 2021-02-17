/*
 *    Copyright (c) 2013-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SMASH_EXPERIMENT_H_
#define SRC_INCLUDE_SMASH_EXPERIMENT_H_

#include <algorithm>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "actionfinderfactory.h"
#include "actions.h"
#include "bremsstrahlungaction.h"
#include "chrono.h"
#include "decayactionsfinder.h"
#include "decayactionsfinderdilepton.h"
#include "energymomentumtensor.h"
#include "fourvector.h"
#include "grandcan_thermalizer.h"
#include "grid.h"
#include "hypersurfacecrossingaction.h"
#include "outputparameters.h"
#include "pauliblocking.h"
#include "potential_globals.h"
#include "potentials.h"
#include "propagation.h"
#include "quantumnumbers.h"
#include "scatteractionphoton.h"
#include "scatteractionsfinder.h"
#include "stringprocess.h"
#include "thermalizationaction.h"
// Output
#include "binaryoutput.h"
#ifdef SMASH_USE_HEPMC
#include "hepmcoutput.h"
#endif
#include "icoutput.h"
#include "oscaroutput.h"
#include "thermodynamicoutput.h"
#ifdef SMASH_USE_ROOT
#include "rootoutput.h"
#endif
#include "vtkoutput.h"
#include "wallcrossingaction.h"

namespace std {
/**
 * Print time span in a human readable way:
 * time < 10 min => seconds
 * 10 min < time < 3 h => minutes
 * time > 3h => hours
 *
 * \note This operator has to be in the \c std namespace for argument dependent
 * lookup to find it. If it were in the smash namespace then the code would not
 * compile since none of its arguments is a type from the smash namespace.
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

namespace smash {
static constexpr int LMain = LogArea::Main::id;
static constexpr int LInitialConditions = LogArea::InitialConditions::id;

/**
 * Non-template interface to Experiment<Modus>.
 *
 * This class allows to call into the public interface of Experiment<Modus>
 * without the need to know the specific `Modus`. The interface is meant for
 * `main()` to set up the experiment and then run takes over.
 */
class ExperimentBase {
 public:
  ExperimentBase() = default;
  /**
   * The virtual destructor avoids undefined behavior when destroying derived
   * objects.
   */
  virtual ~ExperimentBase() = default;

  /**
   * Factory method that creates and initializes a new Experiment<Modus>.
   *
   * This function creates a new Experiment object. The Modus template
   * argument is determined by the \p config argument.
   *
   * \param[in] config The configuration object that sets all initial conditions
   *            of the experiment.
   * \param[in] output_path The directory where the output files are written.
   *
   * \return An owning pointer to the Experiment object, using the
   *         ExperimentBase interface.
   *
   * \throws InvalidModusRequest This exception is thrown if the \p
   *         Modus string in the \p config object does not contain a valid
   *         string.
   *
   * Most of the Configuration values are read starting from this function. The
   * configuration itself is documented in \subpage input_general_
   */
  static std::unique_ptr<ExperimentBase> create(Configuration config,
                                                const bf::path &output_path);

  /**
   * Runs the experiment.
   *
   * The constructor does the setup of the experiment. The run function executes
   * the complete experiment.
   */
  virtual void run() = 0;

  /**
   * \ingroup exception
   * Exception class that is thrown if an invalid modus is requested from the
   * Experiment factory.
   */
  struct InvalidModusRequest : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };
};

template <typename Modus>
class Experiment;
template <typename Modus>
std::ostream &operator<<(std::ostream &out, const Experiment<Modus> &e);

/**
 * The main class, where the simulation of an experiment is executed.
 *
 * The Experiment class owns all data (maybe indirectly) relevant for the
 * execution of the experiment simulation. The experiment can be conducted in
 * different running modi. Since the abstraction of these differences should not
 * incur any overhead, the design is built around the Policy pattern.
 *
 * The Policy pattern was defined by Andrei Alexandrescu in his book "Modern C++
 * Design: Generic Programming and Design Patterns Applied". Addison-Wesley:
 * > A policy defines a class interface or a class template interface.
 * > The interface consists of one or all of the following: inner type
 * > definitions, member functions, and member variables.
 * The policy pattern can also be understood as a compile-time variant of the
 * strategy pattern.
 *
 * The \p Modus template parameter defines the "policy" of the Experiment class.
 * It determines several aspects of the experiment execution *at compile time*.
 * The original strategy pattern would select these differences *at run time*,
 * thus incurring an overhead. This overhead becomes severe in cases where calls
 * to strategy/policy functions are done very frequently. Using the policy
 * pattern, the compiler can fully optimize: It creates a new instance of all
 * functions in Experiment for all different Modus types.
 */
template <typename Modus>
class Experiment : public ExperimentBase {
  friend class ExperimentBase;

 public:
  /**
   * Runs the experiment.
   *
   * The constructor does the setup of the experiment. The run function executes
   * the complete experiment.
   */
  void run() override;

  /**
   * Create a new Experiment.
   *
   * This constructor is only called from the ExperimentBase::create factory
   * method.
   *
   * \param[in] config  The Configuration object contains all initial setup of
   * the experiment. It is forwarded to the constructors of member variables as
   * needed. Note that the object is passed by non-const reference. This is only
   * necessary for bookkeeping: Values are not only read, but actually taken out
   * of the object. Thus, all values that remain were not used. \param[in]
   * output_path The directory where the output files are written.
   */
  explicit Experiment(Configuration config, const bf::path &output_path);

  /**
   * This is called in the beginning of each event. It initializes particles
   * according to selected modus, resets the clock and saves the initial
   * conserved quantities for subsequent sanity checks.
   */
  void initialize_new_event(int event_number);

  /**
   * Runs the time evolution of an event with fixed-sized time steps or without
   * timesteps, from action to actions.
   * Within one timestep (fixed) evolution from action to action
   * is invoked.
   */
  void run_time_evolution();

  /// Performs the final decays of an event
  void do_final_decays();

  /**
   * Output at the end of an event
   *
   * \param[in] evt_num Number of the event
   */
  void final_output(const int evt_num);

  /**
   * Provides external access to SMASH particles. This is helpful if SMASH
   * is used as a 3rd-party library.
   */
  Particles *particles() { return &particles_; }

  /**
   * Provides external access to SMASH calculation modus. This is helpful if
   * SMASH is used as a 3rd-party library.
   */
  Modus *modus() { return &modus_; }

 private:
  /**
   * Perform the given action.
   *
   * \tparam Container type that holds the particles before the action.
   * \param[in] action The action to perform. If it performs, it'll modify
   *                   the private member particles_.
   * \param[in] particles_before_actions A container with the ParticleData
   *                 from this time step before any actions were performed.
   * \return False if the action is rejected either due to invalidity or
   *         Pauli-blocking, or true if it's accepted and performed.
   */
  template <typename Container>
  bool perform_action(Action &action,
                      const Container &particles_before_actions);
  /**
   * Create a list of output files
   *
   * \param[in] format Format of the output file (e.g. Root, Oscar, Vtk)
   * \param[in] content Content of the output (e.g. particles, collisions)
   * \param[in] output_path Path of the output file
   * \param[in] par Output options.(e.g. Extended)
   */
  void create_output(const std::string &format, const std::string &content,
                     const bf::path &output_path, const OutputParameters &par);

  /**
   * Propagate all particles until time to_time without any interactions
   * and shine dileptons.
   *
   * \param[in] to_time Time at the end of propagation [fm/c]
   */
  void propagate_and_shine(double to_time);

  /**
   * Performs all the propagations and actions during a certain time interval
   * neglecting the influence of the potentials. This function is called in
   * either the time stepless cases or the cases with time steps. In a time
   * stepless case, the time interval should be equal to the whole evolution
   * time, while in the case with time step, the intervals are given by the
   * time steps.
   *
   * \param[in, out] actions Actions occur during a certain time interval.
   *                 They provide the ending times of the propagations and
   *                 are updated during the time interval.
   */
  void run_time_evolution_timestepless(Actions &actions);

  /// Intermediate output during an event
  void intermediate_output();

  /// Recompute potentials on lattices if necessary.
  void update_potentials();

  /**
   * Calculate the minimal size for the grid cells such that the
   * ScatterActionsFinder will find all collisions within the maximal
   * transverse distance (which is determined by the maximal cross section).
   *
   * \param[in] dt The current time step size [fm/c]
   * \return The minimal required size of cells
   */
  double compute_min_cell_length(double dt) const {
    return std::sqrt(4 * dt * dt + max_transverse_distance_sqr_);
  }

  /// Shortcut for next output time
  double next_output_time() const {
    return parameters_.outputclock->next_time();
  }

  /**
   * Struct of several member variables.
   * These variables are combined into a struct for efficient input to functions
   * outside of this class.
   */
  ExperimentParameters parameters_;

  /// Structure to precalculate and hold parameters for density computations
  DensityParameters density_param_;

  /**
   * Instance of the Modus template parameter. May store modus-specific data
   * and contains modus-specific function implementations.
   */
  Modus modus_;

  /// Complete particle list
  Particles particles_;

  /**
   * An instance of potentials class, that stores parameters of potentials,
   * calculates them and their gradients.
   */
  std::unique_ptr<Potentials> potentials_;

  /**
   * An instance of PauliBlocker class that stores parameters needed
   * for Pauli blocking calculations and computes phase-space density.
   */
  std::unique_ptr<PauliBlocker> pauli_blocker_;

  /**
   * A list of output formaters. They will be called to write the state of the
   * particles to file.
   */
  OutputsList outputs_;

  /// The Dilepton output
  OutputPtr dilepton_output_;

  /// The Photon output
  OutputPtr photon_output_;

  /**
   * nucleon_has_interacted_ labels whether the particles in the nuclei
   * have experienced any collisions or not. It's only valid in
   * the ColliderModus, so is set as an empty vector by default.
   */
  std::vector<bool> nucleon_has_interacted_ = {};
  /**
   * Whether the projectile and the target collided.
   */
  bool projectile_target_interact_ = false;

  /**
   * The initial nucleons in the ColliderModus propagate with
   * beam_momentum_, if Fermi motion is frozen. It's only valid in
   * the ColliderModus, so is set as an empty vector by default.
   */
  std::vector<FourVector> beam_momentum_ = {};

  /// The Action finder objects
  std::vector<std::unique_ptr<ActionFinderInterface>> action_finders_;

  /// The Dilepton Action Finder
  std::unique_ptr<DecayActionsFinderDilepton> dilepton_finder_;

  /// The (Scatter) Actions Finder for Direct Photons
  std::unique_ptr<ActionFinderInterface> photon_finder_;

  /// Number of fractional photons produced per single reaction
  int n_fractional_photons_;

  /// Baryon density on the lattices
  std::unique_ptr<DensityLattice> jmu_B_lat_;

  /// Isospin projection density on the lattices
  std::unique_ptr<DensityLattice> jmu_I3_lat_;

  /**
   * Custom density on the lattices.
   * In the config user asks for some kind of density for printout.
   * Baryon and isospin projection density are anyway needed for potentials.
   * If user asks for some other density type for printout, it will be handled
   * using jmu_custom variable.
   */
  std::unique_ptr<DensityLattice> jmu_custom_lat_;

  /// Type of density for lattice printout
  DensityType dens_type_lattice_printout_ = DensityType::None;

  /**
   * Lattices for Skyrme potentials (evaluated in the local rest frame) times
   * the baryon flow 4-velocity
   */
  std::unique_ptr<RectangularLattice<FourVector>> UB_lat_ = nullptr;

  /**
   * Lattices for symmetry potentials (evaluated in the local rest frame) times
   * the isospin flow 4-velocity
   */
  std::unique_ptr<RectangularLattice<FourVector>> UI3_lat_ = nullptr;

  /// Lattices for the electric and magnetic components of the Skyrme force
  std::unique_ptr<RectangularLattice<std::pair<ThreeVector, ThreeVector>>>
      FB_lat_;

  /// Lattices for the electric and magnetic component of the symmetry force
  std::unique_ptr<RectangularLattice<std::pair<ThreeVector, ThreeVector>>>
      FI3_lat_;

  /// Lattices of energy-momentum tensors for printout
  std::unique_ptr<RectangularLattice<EnergyMomentumTensor>> Tmn_;

  /// Whether to print the energy-momentum tensor
  bool printout_tmn_ = false;

  /// Whether to print the energy-momentum tensor in Landau frame
  bool printout_tmn_landau_ = false;

  /// Whether to print the 4-velocity in Landau fram
  bool printout_v_landau_ = false;

  /// Whether to print the thermodynamics quantities evaluated on the lattices
  bool printout_lattice_td_ = false;

  /// Instance of class used for forced thermalization
  std::unique_ptr<GrandCanThermalizer> thermalizer_;

  /**
   * Pointer to the string process class object,
   * which is used to set the random seed for PYTHIA objects in each event.
   */
  StringProcess *process_string_ptr_;

  /**
   * Number of events.
   *
   * Event is a single simulation of a physical phenomenon:
   * elementary particle or nucleus-nucleus collision. Result
   * of a single SMASH event is random (by construction)
   * as well as result of one collision in nature. To compare
   * simulation with experiment one has to take ensemble averages,
   * i.e. perform simulation and real experiment many times
   * and compare average results.
   *
   * nevents_ is number of times single phenomenon (particle
   * or nucleus-nucleus collision) will be simulated.
   */
  const int nevents_;

  /// simulation time at which the evolution is stopped.
  const double end_time_;

  /**
   * The clock's timestep size at start up
   *
   * Stored here so that the next event will remember this.
   */
  const double delta_time_startup_;

  /**
   * This indicates whether we force all resonances to decay in the last
   * timestep.
   */
  const bool force_decays_;

  /// This indicates whether to use the grid.
  const bool use_grid_;

  /// This struct contains information on the metric to be used
  const ExpansionProperties metric_;

  /// This indicates whether dileptons are switched on.
  const bool dileptons_switch_;

  /// This indicates whether photons are switched on.
  const bool photons_switch_;

  /// This indicates whether bremsstrahlung is switched on.
  const bool bremsstrahlung_switch_;

  /// This indicates whether the IC output is enabled.
  const bool IC_output_switch_;

  /// This indicates whether to use time steps.
  const TimeStepMode time_step_mode_;

  /// Maximal distance at which particles can interact, squared
  double max_transverse_distance_sqr_ = std::numeric_limits<double>::max();

  /**
   * The conserved quantities of the system.
   *
   * This struct carries the sums of the single particle's various
   * quantities as measured at the beginning of the evolution and can be
   * used to regularly check if they are still good.
   */
  QuantumNumbers conserved_initial_;

  /**
   * The initial total mean field energy in the system.
   * Note: will only be calculated if lattice is on.
   */
  double initial_mean_field_energy_;

  /// system starting time of the simulation
  SystemTimePoint time_start_ = SystemClock::now();

  /// Type of density to be written to collision headers
  DensityType dens_type_ = DensityType::None;

  /**
   *  Total number of interactions for current timestep.
   *  For timestepless mode the whole run time is considered as one timestep.
   */
  uint64_t interactions_total_ = 0;

  /**
   *  Total number of interactions for previous timestep.
   *  For timestepless mode the whole run time is considered as one timestep.
   */
  uint64_t previous_interactions_total_ = 0;

  /**
   *  Total number of wall-crossings for current timestep.
   *  For timestepless mode the whole run time is considered as one timestep.
   */
  uint64_t wall_actions_total_ = 0;

  /**
   *  Total number of wall-crossings for previous timestep.
   *  For timestepless mode the whole run time is considered as one timestep.
   */
  uint64_t previous_wall_actions_total_ = 0;

  /**
   *  Total number of Pauli-blockings for current timestep.
   *  For timestepless mode the whole run time is considered as one timestep.
   */
  uint64_t total_pauli_blocked_ = 0;

  /**
   *  Total number of particles removed from the evolution in
   *  hypersurface crossing actions.
   */
  uint64_t total_hypersurface_crossing_actions_ = 0;

  /**
   *  Total number of discarded interactions, because they were invalidated
   *  before they could be performed.
   */
  uint64_t discarded_interactions_total_ = 0;

  /**
   * Total energy removed from the system in hypersurface crossing actions.
   *
   */
  double total_energy_removed_ = 0.0;

  /// random seed for the next event.
  int64_t seed_ = -1;

  /**
   * \ingroup logging
   * Writes the initial state for the Experiment to the output stream.
   * It automatically appends the output of the current Modus.
   */
  friend std::ostream &operator<<<>(std::ostream &out, const Experiment &e);
};

/// Creates a verbose textual description of the setup of the Experiment.
template <typename Modus>
std::ostream &operator<<(std::ostream &out, const Experiment<Modus> &e) {
  out << "End time: " << e.end_time_ << " fm/c\n";
  out << e.modus_;
  return out;
}

template <typename Modus>
void Experiment<Modus>::create_output(const std::string &format,
                                      const std::string &content,
                                      const bf::path &output_path,
                                      const OutputParameters &out_par) {
  logg[LExperiment].info() << "Adding output " << content << " of format "
                           << format << std::endl;

  if (format == "VTK" && content == "Particles") {
    outputs_.emplace_back(
        make_unique<VtkOutput>(output_path, content, out_par));
  } else if (format == "Root") {
#ifdef SMASH_USE_ROOT
    if (content == "Initial_Conditions") {
      outputs_.emplace_back(
          make_unique<RootOutput>(output_path, "SMASH_IC", out_par));
    } else {
      outputs_.emplace_back(
          make_unique<RootOutput>(output_path, content, out_par));
    }
#else
    logg[LExperiment].error(
        "Root output requested, but Root support not compiled in");
#endif
  } else if (format == "Binary") {
    if (content == "Collisions" || content == "Dileptons" ||
        content == "Photons") {
      outputs_.emplace_back(
          make_unique<BinaryOutputCollisions>(output_path, content, out_par));
    } else if (content == "Particles") {
      outputs_.emplace_back(
          make_unique<BinaryOutputParticles>(output_path, content, out_par));
    } else if (content == "Initial_Conditions") {
      outputs_.emplace_back(make_unique<BinaryOutputInitialConditions>(
          output_path, content, out_par));
    }
  } else if (format == "Oscar1999" || format == "Oscar2013") {
    outputs_.emplace_back(
        create_oscar_output(format, content, output_path, out_par));
  } else if (content == "Thermodynamics" && format == "ASCII") {
    outputs_.emplace_back(
        make_unique<ThermodynamicOutput>(output_path, content, out_par));
  } else if (content == "Thermodynamics" && format == "VTK") {
    printout_lattice_td_ = true;
    outputs_.emplace_back(
        make_unique<VtkOutput>(output_path, content, out_par));
  } else if (content == "Initial_Conditions" && format == "ASCII") {
    outputs_.emplace_back(
        make_unique<ICOutput>(output_path, "SMASH_IC", out_par));
  } else if (content == "HepMC" && format == "ASCII") {
#ifdef SMASH_USE_HEPMC
    outputs_.emplace_back(make_unique<HepMcOutput>(
        output_path, "SMASH_HepMC", out_par, modus_.total_N_number(),
        modus_.proj_N_number()));
#else
    logg[LExperiment].error(
        "HepMC output requested, but HepMC support not compiled in");
#endif
  } else {
    logg[LExperiment].error()
        << "Unknown combination of format (" << format << ") and content ("
        << content << "). Fix the config.";
  }
}

/**
 * Gathers all general Experiment parameters.
 *
 * \param[in, out] config Configuration element
 * \return The ExperimentParameters struct filled with values from the
 *         Configuration
 */
ExperimentParameters create_experiment_parameters(Configuration config);

/*!\Userguide
 * \page input_general_
 * \key End_Time (double, required): \n
 * The time in fm after which the evolution is stopped. Note
 * that the starting time depends on the chosen Modus.
 *
 * \key Randomseed (int, required): \n
 * Initial seed for the random number generator. If this is
 * negative, the seed will be randomly generated by the operating system.
 *
 * \key Nevents (int, required): \n
 * Number of events to calculate.
 *
 * \key Use_Grid (bool, optional, default = true): \n
 * \li \key true - A grid is used to reduce the combinatorics of interaction
 * lookup \n \li \key false - No grid is used.
 *
 * \key Time_Step_Mode (string, optional, default = Fixed): \n
 * The mode of time stepping. Possible values: \n
 * \li \key None - Delta_Time is set to the End_Time.  Cannot be used with
 * potentials. \n
 * \li \key Fixed - Fixed-sized time steps at which collision-finding grid is
 * created.  More efficient for systems with many particles. The Delta_Time is
 * provided by user.\n
 *
 * For Delta_Time explanation see \ref input_general_.
 *
 * \key Metric_Type (string, optional, default = NoExpansion): \n
 * Select which kind of expansion the metric should have. This needs only be
 * specified for the sphere modus:
 * \li \key NoExpansion - Default SMASH run, with Minkowski metric \n
 * \li \key MasslessFRW - FRW expansion going as t^(1/2)
 * \li \key MassiveFRW - FRW expansion going as t^(2/3)
 * \li \key Exponential - FRW expansion going as e^(t/2)
 *
 * \key Expansion_Rate (double, optional, default = 0.1): \n
 * Corresponds to the speed of expansion of the universe in non minkowski
 * metrics if MetricType is any other than \key NoExpansion. \n
 * It corresponds to \f$b_r/l_0\f$ if the metric type is \key MasslessFRW or
 * \key MassiveFRW, and to the parameter b in the Exponential expansion where
 * \f$a(t) ~ e^{bt/2}\f$. \n
 *
 * \page input_collision_term_ Collision Term
 *
 * \key Two_to_One (bool, optional, default = \key true) \n
 * Enable 2 <--> 1 processes (resonance formation and decays).
 *
 * \key Included_2to2 (list of 2 <--> 2 reactions,
 * optional, default = ["All"]) \n
 * List that contains all possible 2 <--> 2 process categories. Each process of
 * the listed category can be performed within the simulation. Possible
 * categories are:
 * \li \key "Elastic" - elastic binary scatterings
 * \li \key "NN_to_NR" - nucleon + nucleon <--> nucleon + resonance
 * \li \key "NN_to_DR" - nucleon + nucleon <--> delta + resonance
 * \li \key "KN_to_KN" - kaon + nucleon <--> kaon + nucleon
 * \li \key "KN_to_KDelta" - kaon + nucleon <--> kaon + delta
 * \li \key "Strangeness_exchange" - processes with strangeness exchange
 * \li \key "NNbar" - annihilation processes, when NNbar_treatment is set to
 *  resonances; this is superseded if NNbar_treatment is set to anything else
 * \li \key "PiDeuteron_to_NN" - deuteron + pion <--> nucleon + nucleon and
 *  its CPT-conjugate
 * \li \key "PiDeuteron_to_pidprime" - deuteron + pion <--> d' + pion
 * \li \key "NDeuteron_to_Ndprime" - deuteron + (anti-)nucleon <-->
 *   d' + (anti-)nucleon, and their CPT-conjugates
 * \li \key "All" - include all binary processes, no necessity to list each
 * single category
 *
 * Detailed balance is preserved by these reaction switches: if a forward
 * reaction is off then the reverse is automatically off too.
 *
 * \key Multi_Particle_Reactions (list of reactions with more than 2 in- or
 * outgoing particles, optional, default = []) \n
 * List that contains all possible multi-particle process categories. Multi
 * particle reactions only work with the stochastic collison criterion.
 * See also example below. Possible categories are:
 * \li \key "Meson_3to1" - Mesonic 3-to-1 reactions:
 * \f$\pi^0\pi^+\pi^-\leftrightarrow\omega\f$,
 * \f$\pi^0\pi^+\pi^-\leftrightarrow\phi\f$,
 * \f$\eta\pi^+\pi^-\leftrightarrow\eta'\f$ and
 * \f$\eta\pi^0\pi^0\leftrightarrow\eta'\f$. Since detailed balance is enforced,
 * the corresponding decays also have to be added in decaymodes.txt to enable
 * the reactions.
 * \li \key "Deuteron_3to2" - Deuteron 3-to-2 reactions:
 * \f$\pi pn\leftrightarrow\pi d\f$, \f$Npn\leftrightarrow Nd\f$ and
 * \f$\bar{N}pn\leftrightarrow \bar{N}d\f$. The deuteron has to be uncommented
 * in particles.txt, too. 2-body reactions involving the d' have to be exluded
 * (no \key "PiDeuteron_to_pidprime" and \key "NDeuteron_to_Ndprime" in
 * \key Included_2to2), since they effectively yield the same reaction.
 * (The same cross section is used as for the d' reactions, therefore the d' in
 * particles.txt and its decay in decaymodest.txt also need to be uncommented.)
 *
 * \key Force_Decays_At_End (bool, optional, default = \key true): \n
 * \li \key true - Force all resonances to decay after last timestep \n
 * \li \key false - Don't force decays (final output can contain resonances)
 *
 * \key No_Collisions (bool, optional, default = false) \n
 * Disable all possible collisions, only allow decays to occur
 * if not forbidden by other options. Useful for running SMASH
 * as a decay afterburner, but not recommended in general, because
 * it breaks the detailed balance.
 *
 * \key NNbar_Treatment (string, optional, default = "strings")
 * \li \key "no annihilation" - No annihilation of NNbar is performed.
 * \li \key "resonances" - Annhilation through NNbar → ρh₁(1170); combined with
 *  ρ → ππ and h₁(1170) → πρ, which gives 5 pions on average. This option
 *  requires "NNbar" to be enabled in Included_2to2.
 * \li \key "strings" - Annihilation through string fragmentation.
 *
 * \key Use_AQM (bool, optional, default = \key true) \n
 * Turn on AQM cross-sections for exotic combination of particles
 * (baryon-baryon cross-sections are scaled from proton-proton high energy
 * parametrization, for example). This includes both elastic and non-elastic
 * contributions; non-elastic contributions go through string fragmentation.
 * Turning off strings or elastic collisions while leaving this on will
 * result in the corresponding part of the AQM cross-sections to also be off.
 * Cross-sections parametrization are scaled according to
 * \f[
 * \frac{\sigma^{AQM}_{\mathrm{process}}}{\sigma^{AQM}_\mathrm{ref\_process}}
 * \sigma^{param}_\mathrm{ref\_process}\f]
 * where \f$ \sigma^{AQM}_x = 40 \left( \frac{2}{3} \right)^{n_{meson}}
 * (1 - 0.4 x^s_1) (1 - 0.4 x^s_2) \f$, with \f$n_{meson}\f$ being the number
 * of mesons in the process, \f$x^s_{1,2}\f$ the fraction of strange quarks in
 * the participant. "process" is then a generic process and "ref_process" a
 * reference process such as PP for which solid parametrizations exist.
 * (\iref{Bass:1998ca})
 *
 * \key Resonance_Lifetime_Modifier (double, optional, default = 1.0)\n
 * Multiplicative factor by which to scale the resonance lifetimes up or down.
 * This additionally has the effect of modifying the initial densities by
 * the same factor in the case of a box initialized with thermal multiplicities
 * (see \key Use_Thermal_Multiplicities). WARNING: This option is not fully
 * physically consistent with some of the other asssumptions used in SMASH;
 * notably, modifying this value WILL break detailed balance in any gas
 * which allows resonances to collide inelastically, as this option breaks the
 * relationship between the width and lifetime of resonances. Note as well that
 * in such gases, using a value of 0.0 is known to make SMASH hang; it is
 * recommended to use a small non-zero value instead in these cases.

 * \key Strings_with_Probability (bool, optional, default = \key true): \n
 * \li \key true - String processes are triggered according to a probability
 *                 increasing smoothly with the collisional energy from 0 to 1
 *                 in a certain energy window. At energies beyond that window,
 *                 all the inelastic scatterings are via strings, while at the
 *                 energies below that window, all the scatterings are via
 *                 non-string processes. One should be careful that in this
 *                 approach, the scatterings via resoances are also suppressed
 *                 in the intermediate energy region, and vanishes at high
 *                 energies, e.g. pπ→Δ→ΣK can't happen at a collisional energy
 *                 beyond 2.2 GeV in this approach. Therefore, the cross
 *                 sections of the scatterings to the certain final states,
 *                 which might be crucial for the production of the rare
 *                 species, will be reduced at the high energies. \n
 * \li \key false - String processes always happen as long as the collisional
 *                  energy exceeds the threshold value by 0.9 GeV, and the
 *                  parametrized total cross section is larger than the sum of
 *                  cross sections contributed by the non-string processes. The
 *                  string cross section is thus obtained by taking the
 *                  difference between them.
 *
 */
template <typename Modus>
Experiment<Modus>::Experiment(Configuration config, const bf::path &output_path)
    : parameters_(create_experiment_parameters(config)),
      density_param_(DensityParameters(parameters_)),
      modus_(config["Modi"], parameters_),
      particles_(),
      nevents_(config.take({"General", "Nevents"})),
      end_time_(config.take({"General", "End_Time"})),
      delta_time_startup_(parameters_.labclock->timestep_duration()),
      force_decays_(
          config.take({"Collision_Term", "Force_Decays_At_End"}, true)),
      use_grid_(config.take({"General", "Use_Grid"}, true)),
      metric_(
          config.take({"General", "Metric_Type"}, ExpansionMode::NoExpansion),
          config.take({"General", "Expansion_Rate"}, 0.1)),
      dileptons_switch_(
          config.take({"Collision_Term", "Dileptons", "Decays"}, false)),
      photons_switch_(config.take(
          {"Collision_Term", "Photons", "2to2_Scatterings"}, false)),
      bremsstrahlung_switch_(
          config.take({"Collision_Term", "Photons", "Bremsstrahlung"}, false)),
      IC_output_switch_(config.has_value({"Output", "Initial_Conditions"})),
      time_step_mode_(
          config.take({"General", "Time_Step_Mode"}, TimeStepMode::Fixed)) {
  logg[LExperiment].info() << *this;

  if (parameters_.coll_crit == CollisionCriterion::Stochastic &&
      time_step_mode_ != TimeStepMode::Fixed) {
    throw std::invalid_argument(
        "The stochastic criterion can only be employed for fixed time step "
        "mode!");
  }

  // create finders
  if (dileptons_switch_) {
    dilepton_finder_ = make_unique<DecayActionsFinderDilepton>();
  }
  if (photons_switch_ || bremsstrahlung_switch_) {
    n_fractional_photons_ =
        config.take({"Collision_Term", "Photons", "Fractional_Photons"}, 100);
  }
  if (parameters_.two_to_one) {
    if (parameters_.res_lifetime_factor < 0.) {
      throw std::invalid_argument(
          "Resonance lifetime modifier cannot be negative!");
    }
    if (parameters_.res_lifetime_factor < really_small) {
      logg[LExperiment].warn(
          "Resonance lifetime set to zero. Make sure resonances cannot "
          "interact",
          "inelastically (e.g. resonance chains), else SMASH is known to "
          "hang.");
    }
    action_finders_.emplace_back(
        make_unique<DecayActionsFinder>(parameters_.res_lifetime_factor));
  }
  bool no_coll = config.take({"Collision_Term", "No_Collisions"}, false);
  if ((parameters_.two_to_one || parameters_.included_2to2.any() ||
       parameters_.included_multi.any() || parameters_.strings_switch) &&
      !no_coll) {
    auto scat_finder = make_unique<ScatterActionsFinder>(
        config, parameters_, nucleon_has_interacted_, modus_.total_N_number(),
        modus_.proj_N_number());
    max_transverse_distance_sqr_ =
        scat_finder->max_transverse_distance_sqr(parameters_.testparticles);
    process_string_ptr_ = scat_finder->get_process_string_ptr();
    action_finders_.emplace_back(std::move(scat_finder));
  } else {
    max_transverse_distance_sqr_ =
        parameters_.maximum_cross_section / M_PI * fm2_mb;
    process_string_ptr_ = NULL;
  }
  if (modus_.is_box()) {
    action_finders_.emplace_back(
        make_unique<WallCrossActionsFinder>(parameters_.box_length));
  }
  if (IC_output_switch_) {
    if (!modus_.is_collider()) {
      throw std::runtime_error(
          "Initial conditions can only be extracted in collider modus.");
    }
    double proper_time;
    if (config.has_value({"Output", "Initial_Conditions", "Proper_Time"})) {
      // Read in proper time from config
      proper_time =
          config.take({"Output", "Initial_Conditions", "Proper_Time"});
    } else {
      // Default proper time is the passing time of the two nuclei
      double default_proper_time = modus_.nuclei_passing_time();
      double lower_bound =
          config.take({"Output", "Initial_Conditions", "Lower_Bound"}, 0.5);
      if (default_proper_time >= lower_bound) {
        proper_time = default_proper_time;
      } else {
        logg[LInitialConditions].warn()
            << "Nuclei passing time is too short, hypersurface proper time set "
               "to tau = "
            << lower_bound << " fm.";
        proper_time = lower_bound;
      }
    }
    action_finders_.emplace_back(
        make_unique<HyperSurfaceCrossActionsFinder>(proper_time));
  }

  if (config.has_value({"Collision_Term", "Pauli_Blocking"})) {
    logg[LExperiment].info() << "Pauli blocking is ON.";
    pauli_blocker_ = make_unique<PauliBlocker>(
        config["Collision_Term"]["Pauli_Blocking"], parameters_);
  }
  // In collider setup with sqrts >= 200 GeV particles don't form continuously
  ParticleData::formation_power_ =
      config.take({"Collision_Term", "Power_Particle_Formation"},
                  modus_.sqrt_s_NN() >= 200. ? -1. : 1.);

  /*!\Userguide
   * \page input_general_
   *
   * \n
   * **Example: Configuring General Properties**\n
   * The following example provides a possibility for the \key General
   * configuration.
   *
   *\verbatim
   General:
       Modus: Collider
       Delta_Time: 0.1
       Testparticles: 1
       Gaussian_Sigma: 1.0
       Gauss_Cutoff_In_Sigma: 3.0
       End_Time: 100.0
       Randomseed: -1
       Nevents: 20
       Use_Grid: True
       Time_Step_Mode: Fixed
   \endverbatim
   *
   * In the case of an expanding sphere setup, change the \key Modus and provide
   * further information about the expansion.
   *\verbatim
       Modus: Sphere
       MetricType: MasslessFRW
       Expansion_Rate: 0.1
   \endverbatim
   *
   **/

  // create outputs
  logg[LExperiment].trace(SMASH_SOURCE_LOCATION, " create OutputInterface objects");

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
   * \subpage configuring_output_.
   *
   * The list of possible contents follows:
   *
   * - \b Particles  List of particles at regular time intervals in the
   *                 computational frame or (optionally) only at the event end.
   *   - Available formats: \ref format_oscar_particlelist,
   *      \ref format_binary_, \ref format_root, \ref format_vtk
   * - \b Collisions List of interactions: collisions, decays, box wall
   *                 crossings and forced thermalizations. Information about
   *                 incoming, outgoing particles and the interaction itself
   *                 is printed out.
   *   - Available formats: \ref format_oscar_collisions, \ref format_binary_,
   *                 \ref format_root
   * - \b Dileptons  Special dilepton output, see \subpage output_dileptons.
   *   - Available formats: \ref format_oscar_collisions,
   *                   \ref format_binary_ and \ref format_root
   * - \b Photons    Special photon output, see \subpage output_photons.
   *   - Available formats: \ref format_oscar_collisions,
   *                   \ref format_binary_ and \ref format_root.
   * - \b Thermodynamics   This output allows to print out thermodynamic
   *          quantities, see \ref Thermodynamics.
   *    - Available formats: \ref thermodyn_output_user_guide_,
   *      \ref output_vtk_lattice_
   * - \b Initial_Conditions  Special initial conditions output, see
   *                          \subpage input_ic for details
   *   - Available formats: \ref format_oscar_particlelist, \ref
   * IC_output_user_guide_
   * - \b HepMC  List of intial and final particles in HepMC3 event record, see
   *                          \subpage hepmc_output_user_guide_ for details
   *   - Available formats: \ref hepmc_output_user_guide_format_
   *
   * \n
   * \anchor list_of_output_formats
   * Output formats
   * --------------
   *
   * For choosing output formats see
   * \ref configuring_output_.
   * Every output content can be printed out in several formats:
   * - \b "Oscar1999", \b "Oscar2013" - human-readable text output\n
   *   - For "Particles" content: \subpage format_oscar_particlelist
   *   - For "Collisions" content: \subpage format_oscar_collisions
   *   - General block structure of OSCAR formats: \subpage oscar_general_
   * - \b "Binary" - binary, not human-readable output
   *   - Faster to read and write than text outputs
   *   - Saves coordinates and momenta with the full double precision
   *   - General file structure is similar to \ref oscar_general_
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
   *   - Used for "Thermodynamics", "Initial_Conditions" and "HepMC", see
   * \subpage thermodyn_output_user_guide_
   * \subpage IC_output_user_guide_
   * \ref hepmc_output_user_guide_
   *
   * \note Output of coordinates for the "Collisions" content in
   *       the periodic box has a feature:
   *       \subpage collisions_output_in_box_modus_
   */

  /*!\Userguide
   * \page output_dileptons Dileptons
   * The existence of a dilepton subsection in the collision term section of the
   * configuration file enables the dilepton production. In addition, the
   * dilepton output also needs to be enabled in the output section and dilepton
   * decays have to be uncommented in the used decaymodes.txt file. The output
   * file named Dileptons (followed by the appropriate suffix) is generated when
   * SMASH is executed. It's format is identical to the collision output (see
   * \ref format_oscar_collisions), it does however only contain information
   * about the dilepton decays. \n Further, the block headers differ from the
   * usual collision output: <div class="fragment"> <div class="line"> <span
   * class="preprocessor">
   *  \# interaction in nin out nout rho density weight shining_weight partial
   *  part_weight type proc_type </span></div>
   * </div>
   * where \li \key nin: Number of ingoing
   * particles (initial state particles) \li \key nout: Number of outgoing
   * particles (finalstate particles) \li \key density: Density at the
   * interaction point \li \key shining_weight: Shining weight of the
   * interaction. Explanation follows below. \li \key part_weight: The partial
   * weight of the interaction. For the dileptons, this coincides with the
   * branching ratio. \li \key proc_type: The type of the underlying process.
   * See process_type for possible types.
   *
   * Note, that "interaction", "in", "out", "rho", "weight", "partial" and
   * "type" are no variables, but words that are printed. \n
   * The dilepton output is available in binary, OSCAR1999, OSCAR2013 and
   * OSCAR2013 extended format. \n
   *
   * \n
   * \note
   * As dileptons are treated perturbatively, the produced dileptons are
   * only written to the dilepton output, but neither to the usual collision
   * output, nor to the particle lists.
   **/

  /*!\Userguide
   * \page output_photons Photons
   * The existance of a photon subsection in the output section of the
   * configuration file enables the photon output.
   * If photons are enabled, the output file named Photons (followed by the
   * appropriate suffix) is generated when SMASH is executed. It's format is
   * identical to the collision output (see \ref format_oscar_collisions),
   * it does however only contain information about all particles participating
   * in the photon producing interactions at each timestep. \n
   * Further, the block headers differ from the usual collision output:
   * <div class="fragment">
   * <div class="line"> <span class="preprocessor">
   *  \# interaction in nin out nout rho density weight photon_weight partial
   *  part_weight type proc_type </span></div>
   * </div>
   * where
   * \li \key density: Density at the interaction point
   * \li \key photon_weight: Weight of the photon process relative to the
   * underlying hadonic interaction. Make sure to weigh each photon in your
   * analysis with this value. Otherwise the photon production is highly
   * overestimated.
   * \li \key part_weight: Always 0.0 for photon processes, as they
   * are hardcoded.
   * \li \key proc_type: The type of the underlying process. See
   * \ref process_type for possible types.
   *
   * Note, that "interaction", "in", "out", "rho", "weight", "partial" and
   * "type" are no variables, but words that are printed. \n
   * The photon output is available in binary, OSCAR1999, OSCAR2013 and
   * OSCAR2013 extended format. \n
   *
   **/

  /*!\Userguide
   * \page input_ic Initial Conditions
   * The existance of an initial conditions subsection in the output section of
   * the configuration file enables the IC output. In addition, all particles
   * that cross the hypersurface of predefined proper time are removed from the
   * evolution. This proper time is taken from the \key Proper_Time field
   * in the \key Initial_Conditions subsection when configuring the output. If
   * this information
   * is not provided, the default proper time corresponds to the passing time
   * of the two nuclei, where all primary interactions are expected to have
   * occured: \f[ \tau_0 =
   * (r_\mathrm{p} \ + \ r_\mathrm{t}) \ \left(\left(\frac{\sqrt{s_\mathrm{NN}}}
   * {2 \ m_\mathrm{N}}\right)^2
   * - 1\right)^{-1/2} \f]
   * Therein, \f$ r_\mathrm{p} \f$ and \f$ r_\mathrm{t} \f$ denote the radii of
   * the projectile and target nucleus, respectively, \f$
   * \sqrt{s_\mathrm{NN}}\f$
   * is the collision energy per nucleon and \f$ m_\mathrm{N} \f$ the nucleon
   * mass. Note though that, if the passing time is smaller than 0.5 fm, the
   * default porper time of the hypersurface is taken to be \f$\tau = 0.5 \f$
   * as a minimum bound to ensure the proper time is large enough
   * to also extract reasonable initial conditions at RHIC/LHC energies. If
   * desired, this lowest possible value can also be specifie in the
   * configuration file in the \key Lower_Bound field. \n Once
   * initial conditions are enabled, the output file named SMASH_IC (followed by
   * the appropriate suffix) is generated when SMASH is executed. \n The output
   * is available in Oscar1999, Oscar2013, binary and ROOT format, as well as in
   * an aditional ASCII format (see \ref IC_output_user_guide_). The latter is
   * meant to directly serve
   * as an input for the vHLLE hydrodynamics code (I. Karpenko, P. Huovinen, M.
   * Bleicher: Comput. Phys. Commun. 185, 3016 (2014)).\n \n
   * ### Oscar output
   * In case
   * of the Oscar1999 and Oscar2013 format, the structure is identical to the
   * Oscar Particles Format (see \ref format_oscar_particlelist). \n
   * In contrast
   * to the usual particles output however, the initial conditions output
   * provides a
   * **list of all particles removed from the evolution** at the time when
   * crossing the hypersurface. This implies that neither the initial particle
   * list nor the particle list at each time step is printed.\n The general
   * Oscar structure as described in \ref format_oscar_particlelist is
   * preserved. \n
   * \n
   * ### Binary output
   * The binary initial
   * conditions output also provides a list of all particles removed from the
   * evolution at the time when crossing the hypersurface. For each removed
   * particle a 'p' block is created stores the particle data. The general
   * binary output structure as described in \ref format_binary_ is preserved.\n
   * \n
   * ### ROOT output
   * The initial conditions output in shape of a list of all particles removed
   * from the SMASH evolution when crossing the hypersurface is also available
   * in ROOT format. Neither the initial nor the final particle lists are
   * printed, but the general structure for particle TTrees, as described in
   * \ref format_root, is preserved.
   */
  dens_type_ = config.take({"Output", "Density_Type"}, DensityType::None);
  logg[LExperiment].debug()
      << "Density type printed to headers: " << dens_type_;

  const OutputParameters output_parameters(std::move(output_conf));

  std::vector<std::string> output_contents = output_conf.list_upmost_nodes();
  for (const auto &content : output_contents) {
    auto this_output_conf = output_conf[content.c_str()];
    const std::vector<std::string> formats = this_output_conf.take({"Format"});
    if (output_path == "") {
      continue;
    }
    for (const auto &format : formats) {
      create_output(format, content, output_path, output_parameters);
    }
  }

  /* We can take away the Fermi motion flag, because the collider modus is
   * already initialized. We only need it when potentials are enabled, but we
   * always have to take it, otherwise SMASH will complain about unused
   * options. We have to provide a default value for modi other than Collider.
   */
  const FermiMotion motion =
      config.take({"Modi", "Collider", "Fermi_Motion"}, FermiMotion::Off);
  if (config.has_value({"Potentials"})) {
    if (time_step_mode_ == TimeStepMode::None) {
      logg[LExperiment].error() << "Potentials only work with time steps!";
      throw std::invalid_argument("Can't use potentials without time steps!");
    }
    if (motion == FermiMotion::Frozen) {
      logg[LExperiment].error()
          << "Potentials don't work with frozen Fermi momenta! "
             "Use normal Fermi motion instead.";
      throw std::invalid_argument(
          "Can't use potentials "
          "with frozen Fermi momenta!");
    }
    logg[LExperiment].info() << "Potentials are ON. Timestep is "
                             << parameters_.labclock->timestep_duration();
    // potentials need testparticles and gaussian sigma from parameters_
    potentials_ = make_unique<Potentials>(config["Potentials"], parameters_);
  }

  /*!\Userguide
   * \page input_lattice_ Lattice
   * It is possible to configure a Lattice for the 3D space, which can be
   * useful to speed up the computation of the potentials. Note though, that
   * this goes in hand with a loss of accuracy: If the lattice is
   * applied, the evaluation of the potentials is carried out only on the nodes
   * of the lattice. Intermediate values are interpolated. \n
   * The configuration of a lattice is usually not necessary, it is however
   * required if the Thermodynamic VTK Output (see
   * \ref output_vtk_lattice_) or the \key Potentials_Affect_Thresholds option
   * is enabled. \n
   * The following parameters are only required, if the \key Lattice section is
   * used in the configuration file. Otherwise, no lattice will be used at all.
   *
   *
   * \key Sizes (array<double,3>, required, no default): \n
   *      Sizes of lattice in x, y, z directions in fm.
   *
   * \key Cell_Number (array<int,3>, required, no default): \n
   *      Number of cells in x, y, z directions.
   *
   * \key Origin (array<double,3>, required, no default): \n
   *      Coordinates of the left, down, near corner of the lattice in fm.
   *
   * \key Periodic (bool, required, no default): \n
   *      Use periodic continuation or not. With periodic continuation
   *      x + i * lx is equivalent to x, same for y, z.
   *
   * \key Potentials_Affect_Thresholds (bool, optional, default = false): \n
   * Include potential effects, since mean field potentials change the threshold
   * energies of the actions.
   *
   * For information on the format of the lattice output see
   * \ref output_vtk_lattice_. To configure the
   * thermodynamic output, see \ref input_output_options_.
   *
   * \n
   * Examples: Configuring the Lattice
   * --------------
   * The following example configures the lattice with the origin in (0,0,0),
   * 20 cells of 10 fm size in each direction and with periodic boundary
   * conditions. The potential effects on the thresholds are taken into
   * consideration. Note that, as the origin is by definition the left down near
   * corner of the cell, center is located at (5, 5, 5).
   *
   *\verbatim
   Lattice:
       Origin:    [0.0, 0.0, 0.0]
       Sizes:    [10.0, 10.0, 10.0]
       Cell_Number:    [20, 20, 20]
       Periodic: True
       Potentials_Affect_Thresholds: True
   \endverbatim
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
        UB_lat_ = make_unique<RectangularLattice<FourVector>>(
            l, n, origin, periodic, LatticeUpdate::EveryTimestep);
        FB_lat_ = make_unique<
            RectangularLattice<std::pair<ThreeVector, ThreeVector>>>(
            l, n, origin, periodic, LatticeUpdate::EveryTimestep);
      }
      if (potentials_->use_symmetry()) {
        jmu_I3_lat_ = make_unique<DensityLattice>(l, n, origin, periodic,
                                                  LatticeUpdate::EveryTimestep);
        UI3_lat_ = make_unique<RectangularLattice<FourVector>>(
            l, n, origin, periodic, LatticeUpdate::EveryTimestep);
        FI3_lat_ = make_unique<
            RectangularLattice<std::pair<ThreeVector, ThreeVector>>>(
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
    logg[LExperiment].error(
        "If you want Thermodynamic VTK output, configure a lattice for it.");
  }
  // Warning for the mean field calculation if lattice is not on.
  if ((potentials_ != nullptr) && (jmu_B_lat_ == nullptr)) {
    logg[LExperiment].warn() << "Lattice is NOT used. Mean fields are "
                             << "not going to be calculated.";
  }

  // Store pointers to potential and lattice accessible for Action
  if (parameters_.potential_affect_threshold) {
    UB_lat_pointer = UB_lat_.get();
    UI3_lat_pointer = UI3_lat_.get();
    pot_pointer = potentials_.get();
  }

  // Create forced thermalizer
  if (config.has_value({"Forced_Thermalization"})) {
    Configuration &&th_conf = config["Forced_Thermalization"];
    thermalizer_ = modus_.create_grandcan_thermalizer(th_conf);
  }

  /* Take the seed setting only after the configuration was stored to a file
   * in smash.cc */
  seed_ = config.take({"General", "Randomseed"});
}

/// String representing a horizontal line.
const std::string hline(113, '-');

/**
 * Generate the tabulated string which will be printed to the screen when
 * SMASH is running
 *
 * \param[in] particles The interacting particles. Their information will be
 *            used to check the conservation of the total energy and momentum.
 *	      the total number of the particles will be used and printed as
 *	      well.
 * \param[in] scatterings_this_interval Number of the scatterings occur within
 *            the current timestep.
 * \param[in] conserved_initial Initial quantum numbers needed to check the
 *            conservations.
 * \param[in] time_start Moment in the REAL WORLD when SMASH starts to run [s].
 * \param[in] time Current moment in SMASH [fm/c].
 * \param[in] E_mean_field Value of the mean-field contribution to the total
 *            energy of the system at the current time.
 * \param[in] E_mean_field_initial Value of the mean-field contribution to the
 *            total energy of the system at t=0.
 * \return 'Current time in SMASH [fm/c]', 'Total kinetic energy in the system
 *         [GeV]', 'Total mean field energy in the system [GeV]', 'Total energy
 *         in the system [GeV]', 'Total energy per particle [GeV]', 'Deviation
 *         of the energy per particle from the initial value [GeV]', 'Number of
 *         scatterings that occurred within the timestep', 'Total particle
 *         number', 'Computing time consumed'.
 */
std::string format_measurements(const Particles &particles,
                                uint64_t scatterings_this_interval,
                                const QuantumNumbers &conserved_initial,
                                SystemTimePoint time_start, double time,
                                double E_mean_field,
                                double E_mean_field_initial);
/**
 * Calculate the total mean field energy of the system; this will be printed to
 * the screen when SMASH is running. Using the baryon density lattice is
 * necessary.
 *
 * \param[in] potentials Parameters of the potentials used in the simulation.
 * \param[in] jmu_B_lat Lattice of baryon density and baryon current values as
 *            well as their gradients at each lattice node.
 * \param[in] parameters Parameters of the experiment, needed for the access to
 *            the number of testparticles.
 * \return Total mean field energy in the Box.
 */
double calculate_mean_field_energy(
    const Potentials &potentials,
    RectangularLattice<smash::DensityOnLattice> &jmu_B_lat,
    const ExperimentParameters &parameters);

/**
 * Generate the EventInfo object which is passed to outputs_.
 *
 * \param[in] particles The interacting particles. Their information will be
 *            used to check the conservation of the total energy and momentum.
 *	      the total number of the particles will be used and printed as
 *	      well.
 * \param[in] E_mean_field Value of the mean-field contribution to the total
 *            energy of the system at the current time.
 * \param[in] modus_impact_parameter The impact parameter
 * \param[in] parameters structure that holds various global parameters
 *            such as testparticle number, see \ref ExperimentParameters
 * \param[in] projectile_target_interact true if there was at least one
 *            collision
 */
EventInfo fill_event_info(const Particles &particles, double E_mean_field,
                          double modus_impact_parameter,
                          const ExperimentParameters &parameters,
                          bool projectile_target_interact);

template <typename Modus>
void Experiment<Modus>::initialize_new_event(int event_number) {
  random::set_seed(seed_);
  logg[LExperiment].info() << "random number seed: " << seed_;
  /* Set seed for the next event. It has to be positive, so it can be entered
   * in the config.
   *
   * We have to be careful about the minimal integer, whose absolute value
   * cannot be represented. */
  int64_t r = random::advance();
  while (r == INT64_MIN) {
    r = random::advance();
  }
  seed_ = std::abs(r);
  /* Set the random seed used in PYTHIA hadronization
   * to be same with the SMASH one.
   * In this way we ensure that the results are reproducible
   * for every event if one knows SMASH random seed. */
  if (process_string_ptr_ != NULL) {
    process_string_ptr_->init_pythia_hadron_rndm();
  }

  particles_.reset();

  // Sample particles according to the initial conditions
  double start_time = modus_.initial_conditions(&particles_, parameters_);
  /* For box modus make sure that particles are in the box. In principle, after
   * a correct initialization they should be, so this is just playing it safe.
   */
  modus_.impose_boundary_conditions(&particles_, outputs_);
  // Reset the simulation clock
  double timestep = delta_time_startup_;

  switch (time_step_mode_) {
    case TimeStepMode::Fixed:
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
  std::unique_ptr<UniformClock> clock_for_this_event;
  if (modus_.is_list() && (timestep < 0.0)) {
    throw std::runtime_error(
        "Timestep for the given event is negative. \n"
        "This might happen if the formation times of the input particles are "
        "larger than the specified end time of the simulation.");
  }
  clock_for_this_event = make_unique<UniformClock>(start_time, timestep);
  parameters_.labclock = std::move(clock_for_this_event);

  // Reset the output clock
  parameters_.outputclock->reset(start_time, true);
  // remove time before starting time in case of custom output times.
  parameters_.outputclock->remove_times_in_past(start_time);

  logg[LExperiment].debug(
      "Lab clock: t_start = ", parameters_.labclock->current_time(),
      ", dt = ", parameters_.labclock->timestep_duration());

  /* Save the initial conserved quantum numbers and total momentum in
   * the system for conservation checks */
  conserved_initial_ = QuantumNumbers(particles_);
  wall_actions_total_ = 0;
  previous_wall_actions_total_ = 0;
  interactions_total_ = 0;
  previous_interactions_total_ = 0;
  discarded_interactions_total_ = 0;
  total_pauli_blocked_ = 0;
  projectile_target_interact_ = false;
  total_hypersurface_crossing_actions_ = 0;
  total_energy_removed_ = 0.0;
  // Print output headers
  logg[LExperiment].info() << hline;
  logg[LExperiment].info() << "Time[fm]   Ekin[GeV]   E_MF[GeV]  ETotal[GeV]  "
                           << "ETot/N[GeV]  D(ETot/N)[GeV] Scatt&Decays  "
                           << "Particles     Comp.Time";
  logg[LExperiment].info() << hline;
  double E_mean_field = 0.0;
  if (potentials_) {
    update_potentials();
    // using the lattice is necessary
    if ((jmu_B_lat_ != nullptr)) {
      E_mean_field =
          calculate_mean_field_energy(*potentials_, *jmu_B_lat_, parameters_);
    }
  }
  initial_mean_field_energy_ = E_mean_field;
  logg[LExperiment].info() << format_measurements(
      particles_, 0u, conserved_initial_, time_start_,
      parameters_.labclock->current_time(), E_mean_field,
      initial_mean_field_energy_);

  auto event_info =
      fill_event_info(particles_, E_mean_field, modus_.impact_parameter(),
                      parameters_, projectile_target_interact_);

  // Output at event start
  for (const auto &output : outputs_) {
    output->at_eventstart(particles_, event_number, event_info);
  }
}

template <typename Modus>
template <typename Container>
bool Experiment<Modus>::perform_action(
    Action &action, const Container &particles_before_actions) {
  // Make sure to skip invalid and Pauli-blocked actions.
  if (!action.is_valid(particles_)) {
    discarded_interactions_total_++;
    logg[LExperiment].debug(~einhard::DRed(), "✘ ", action,
                            " (discarded: invalid)");
    return false;
  }
  action.generate_final_state();
  logg[LExperiment].debug("Process Type is: ", action.get_type());
  if (pauli_blocker_ && action.is_pauli_blocked(particles_, *pauli_blocker_)) {
    total_pauli_blocked_++;
    return false;
  }
  if (modus_.is_collider()) {
    /* Mark incoming nucleons as interacted - now they are permitted
     * to collide with nucleons from their native nucleus */
    bool incoming_projectile = false;
    bool incoming_target = false;
    for (const auto &incoming : action.incoming_particles()) {
      assert(incoming.id() >= 0);
      if (incoming.id() < modus_.total_N_number()) {
        nucleon_has_interacted_[incoming.id()] = true;
      }
      if (incoming.id() < modus_.proj_N_number()) {
        incoming_projectile = true;
      }
      if (incoming.id() >= modus_.proj_N_number() &&
          incoming.id() < modus_.total_N_number()) {
        incoming_target = true;
      }
    }
    // Check whether particles from different nuclei interacted.
    if (incoming_projectile & incoming_target) {
      projectile_target_interact_ = true;
    }
  }
  /* Make sure to pick a non-zero integer, because 0 is reserved for "no
   * interaction yet". */
  const auto id_process = static_cast<uint32_t>(interactions_total_ + 1);
  action.perform(&particles_, id_process);
  interactions_total_++;
  if (action.get_type() == ProcessType::Wall) {
    wall_actions_total_++;
  }
  if (action.get_type() == ProcessType::HyperSurfaceCrossing) {
    total_hypersurface_crossing_actions_++;
    total_energy_removed_ += action.incoming_particles()[0].momentum().x0();
  }
  // Calculate Eckart rest frame density at the interaction point
  double rho = 0.0;
  if (dens_type_ != DensityType::None) {
    const FourVector r_interaction = action.get_interaction_point();
    constexpr bool compute_grad = false;
    const bool smearing = true;
    rho = std::get<0>(current_eckart(r_interaction.threevec(),
                                     particles_before_actions, density_param_,
                                     dens_type_, compute_grad, smearing));
  }
  /*!\Userguide
   * \page collisions_output_in_box_modus_ Collision Output in Box Modus
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
      if (output->is_IC_output() &&
          action.get_type() == ProcessType::HyperSurfaceCrossing) {
        output->at_interaction(action, rho);
      } else if (!output->is_IC_output()) {
        output->at_interaction(action, rho);
      }
    }
  }

  // At every collision photons can be produced.
  // Note: We rely here on the lazy evaluation of the arguments to if.
  // It may happen that in a wall-crossing-action sqrt_s raises an exception.
  // Therefore we first have to check if the incoming particles can undergo
  // an em-interaction.
  if (photons_switch_ &&
      ScatterActionPhoton::is_photon_reaction(action.incoming_particles()) &&
      ScatterActionPhoton::is_kinematically_possible(
          action.sqrt_s(), action.incoming_particles())) {
    /* Time in the action constructor is relative to
     * current time of incoming */
    constexpr double action_time = 0.;
    ScatterActionPhoton photon_act(action.incoming_particles(), action_time,
                                   n_fractional_photons_,
                                   action.get_total_weight());

    /**
     * Add a completely dummy process to the photon action. The only important
     * thing is that its cross-section is equal to the cross-section of the
     * hadronic action. This can be done, because the photon action is never
     * actually performed, only the final state is generated and printed to
     * the photon output.
     * Note: The cross_section_scaling_factor can be neglected here, since it
     * cancels out for the weighting, where a ratio of (unscaled) photon
     * cross section and (unscaled) hadronic cross section is taken.
     */
    photon_act.add_dummy_hadronic_process(action.get_total_weight());

    // Now add the actual photon reaction channel.
    photon_act.add_single_process();

    photon_act.perform_photons(outputs_);
  }

  if (bremsstrahlung_switch_ &&
      BremsstrahlungAction::is_bremsstrahlung_reaction(
          action.incoming_particles())) {
    /* Time in the action constructor is relative to
     * current time of incoming */
    constexpr double action_time = 0.;

    BremsstrahlungAction brems_act(action.incoming_particles(), action_time,
                                   n_fractional_photons_,
                                   action.get_total_weight());

    /**
     * Add a completely dummy process to the bremsstrahlung action. The only
     * important thing is that its cross-section is equal to the cross-section
     * of the hadronic action. This can be done, because the bremsstrahlung
     * action is never actually performed, only the final state is generated and
     * printed to the photon output. Note: The cross_section_scaling_factor can
     * be neglected here, since it cancels out for the weighting, where a ratio
     * of (unscaled) photon cross section and (unscaled) hadronic cross section
     * is taken.
     */

    brems_act.add_dummy_hadronic_process(action.get_total_weight());

    // Now add the actual bremsstrahlung reaction channel.
    brems_act.add_single_process();

    brems_act.perform_bremsstrahlung(outputs_);
  }

  logg[LExperiment].debug(~einhard::Green(), "✔ ", action);
  return true;
}

template <typename Modus>
void Experiment<Modus>::run_time_evolution() {
  Actions actions;

  while (parameters_.labclock->current_time() < end_time_) {
    const double t = parameters_.labclock->current_time();
    const double dt =
        std::min(parameters_.labclock->timestep_duration(), end_time_ - t);
    logg[LExperiment].debug("Timestepless propagation for next ", dt, " fm/c.");

    // Perform forced thermalization if required
    if (thermalizer_ &&
        thermalizer_->is_time_to_thermalize(parameters_.labclock)) {
      const bool ignore_cells_under_treshold = true;
      thermalizer_->update_thermalizer_lattice(particles_, density_param_,
                                               ignore_cells_under_treshold);
      const double current_t = parameters_.labclock->current_time();
      thermalizer_->thermalize(particles_, current_t,
                               parameters_.testparticles);
      ThermalizationAction th_act(*thermalizer_, current_t);
      if (th_act.any_particles_thermalized()) {
        perform_action(th_act, particles_);
      }
    }

    if (particles_.size() > 0 && action_finders_.size() > 0) {
      /* (1.a) Create grid. */
      double min_cell_length = compute_min_cell_length(dt);
      logg[LExperiment].debug("Creating grid with minimal cell length ",
                              min_cell_length);
      const auto &grid =
          use_grid_ ? modus_.create_grid(particles_, min_cell_length, dt)
                    : modus_.create_grid(particles_, min_cell_length, dt,
                                         CellSizeStrategy::Largest);

      const double gcell_vol = grid.cell_volume();

      /* (1.b) Iterate over cells and find actions. */
      grid.iterate_cells(
          [&](const ParticleList &search_list) {
            for (const auto &finder : action_finders_) {
              actions.insert(finder->find_actions_in_cell(
                  search_list, dt, gcell_vol, beam_momentum_));
            }
          },
          [&](const ParticleList &search_list,
              const ParticleList &neighbors_list) {
            for (const auto &finder : action_finders_) {
              actions.insert(finder->find_actions_with_neighbors(
                  search_list, neighbors_list, dt, beam_momentum_));
            }
          });
    }

    /* \todo (optimizations) Adapt timestep size here */

    /* (2) Propagation from action to action until the end of timestep */
    run_time_evolution_timestepless(actions);

    /* (3) Update potentials (if computed on the lattice) and
     *     compute new momenta according to equations of motion */
    if (potentials_) {
      update_potentials();
      update_momenta(&particles_, parameters_.labclock->timestep_duration(),
                     *potentials_, FB_lat_.get(), FI3_lat_.get());
    }

    /* (4) Expand universe if non-minkowskian metric; updates
     *     positions and momenta according to the selected expansion */
    if (metric_.mode_ != ExpansionMode::NoExpansion) {
      expand_space_time(&particles_, parameters_, metric_);
    }

    ++(*parameters_.labclock);

    /* (5) Check conservation laws.
     *
     * Check conservation of conserved quantities if potentials and string
     * fragmentation are off.  If potentials are on then momentum is conserved
     * only in average.  If string fragmentation is on, then energy and
     * momentum are only very roughly conserved in high-energy collisions. */
    if (!potentials_ && !parameters_.strings_switch &&
        metric_.mode_ == ExpansionMode::NoExpansion && !IC_output_switch_) {
      std::string err_msg = conserved_initial_.report_deviations(particles_);
      if (!err_msg.empty()) {
        logg[LExperiment].error() << err_msg;
        throw std::runtime_error("Violation of conserved quantities!");
      }
    }
  }

  if (pauli_blocker_) {
    logg[LExperiment].info(
        "Interactions: Pauli-blocked/performed = ", total_pauli_blocked_, "/",
        interactions_total_ - wall_actions_total_);
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

/**
 * Make sure `interactions_total` can be represented as a 32-bit integer.
 * This is necessary for converting to a `id_process`. The latter is 32-bit
 * integer, because it is written like this to binary output.
 *
 * \param[in] interactions_total Total interaction number
 */
inline void check_interactions_total(uint64_t interactions_total) {
  constexpr uint64_t max_uint32 = std::numeric_limits<uint32_t>::max();
  if (interactions_total >= max_uint32) {
    throw std::runtime_error("Integer overflow in total interaction number!");
  }
}

template <typename Modus>
void Experiment<Modus>::run_time_evolution_timestepless(Actions &actions) {
  const double start_time = parameters_.labclock->current_time();
  const double end_time =
      std::min(parameters_.labclock->next_time(), end_time_);
  double time_left = end_time - start_time;
  logg[LExperiment].debug(
      "Timestepless propagation: ", "Actions size = ", actions.size(),
      ", start time = ", start_time, ", end time = ", end_time);

  // iterate over all actions
  while (!actions.is_empty()) {
    // get next action
    ActionPtr act = actions.pop();
    if (!act->is_valid(particles_)) {
      discarded_interactions_total_++;
      logg[LExperiment].debug(~einhard::DRed(), "✘ ", act,
                              " (discarded: invalid)");
      continue;
    }
    if (act->time_of_execution() > end_time) {
      logg[LExperiment].error(
          act, " scheduled later than end time: t_action[fm/c] = ",
          act->time_of_execution(), ", t_end[fm/c] = ", end_time);
    }
    logg[LExperiment].debug(~einhard::Green(), "✔ ", act);

    while (next_output_time() <= act->time_of_execution()) {
      logg[LExperiment].debug("Propagating until output time: ",
                              next_output_time());
      propagate_and_shine(next_output_time());
      ++(*parameters_.outputclock);
      intermediate_output();
    }

    /* (1) Propagate to the next action. */
    logg[LExperiment].debug("Propagating until next action ", act,
                            ", action time = ", act->time_of_execution());
    propagate_and_shine(act->time_of_execution());

    /* (2) Perform action.
     *
     * Update the positions of the incoming particles, because the information
     * in the action object will be outdated as the particles have been
     * propagated since the construction of the action. */
    act->update_incoming(particles_);
    const bool performed = perform_action(*act, particles_);

    /* No need to update actions for outgoing particles
     * if the action is not performed. */
    if (!performed) {
      continue;
    }

    /* (3) Update actions for newly-produced particles. */

    time_left = end_time - act->time_of_execution();
    const ParticleList &outgoing_particles = act->outgoing_particles();
    // Grid cell volume set to zero, since there is no grid
    const double gcell_vol = 0.0;
    for (const auto &finder : action_finders_) {
      // Outgoing particles can still decay, cross walls...
      actions.insert(finder->find_actions_in_cell(outgoing_particles, time_left,
                                                  gcell_vol, beam_momentum_));
      // ... and collide with other particles.
      actions.insert(finder->find_actions_with_surrounding_particles(
          outgoing_particles, particles_, time_left, beam_momentum_));
    }

    check_interactions_total(interactions_total_);
  }

  while (next_output_time() <= end_time) {
    logg[LExperiment].debug("Propagating until output time: ",
                            next_output_time());
    propagate_and_shine(next_output_time());
    ++(*parameters_.outputclock);
    // Avoid duplicating printout at event end time
    if (parameters_.outputclock->current_time() < end_time_) {
      intermediate_output();
    }
  }
  logg[LExperiment].debug("Propagating to time ", end_time);
  propagate_and_shine(end_time);
}

template <typename Modus>
void Experiment<Modus>::intermediate_output() {
  const uint64_t wall_actions_this_interval =
      wall_actions_total_ - previous_wall_actions_total_;
  previous_wall_actions_total_ = wall_actions_total_;
  const uint64_t interactions_this_interval = interactions_total_ -
                                              previous_interactions_total_ -
                                              wall_actions_this_interval;
  previous_interactions_total_ = interactions_total_;
  double E_mean_field = 0.0;
  if (potentials_) {
    // using the lattice is necessary
    if ((jmu_B_lat_ != nullptr)) {
      E_mean_field =
          calculate_mean_field_energy(*potentials_, *jmu_B_lat_, parameters_);
      /*
       * Mean field calculated in a box should remain approximately constant if
       * the system is in equilibrium, and so deviations from its original value
       * may signal a phase transition or other dynamical process. This
       * comparison only makes sense in the Box Modus, hence the condition.
       */
      if (modus_.is_box()) {
        double tmp = (E_mean_field - initial_mean_field_energy_) /
                     (E_mean_field + initial_mean_field_energy_);
        /*
         * This is displayed when the system evolves away from its initial
         * configuration (which is when the total mean field energy in the box
         *  deviates from its initial value).
         */
        if (std::abs(tmp) > 0.01) {
          logg[LExperiment].info()
              << "\n\n\n\t The mean field at t = "
              << parameters_.outputclock->current_time()
              << " [fm/c] differs from the mean field at t = 0:"
              << "\n\t\t                 initial_mean_field_energy_ = "
              << initial_mean_field_energy_ << " [GeV]"
              << "\n\t\t abs[(E_MF - E_MF(t=0))/(E_MF + E_MF(t=0))] = "
              << std::abs(tmp)
              << "\n\t\t                             E_MF/E_MF(t=0) = "
              << E_mean_field / initial_mean_field_energy_ << "\n\n";
        }
      }
    }
  }

  logg[LExperiment].info() << format_measurements(
      particles_, interactions_this_interval, conserved_initial_, time_start_,
      parameters_.outputclock->current_time(), E_mean_field,
      initial_mean_field_energy_);
  const LatticeUpdate lat_upd = LatticeUpdate::AtOutput;

  auto event_info =
      fill_event_info(particles_, E_mean_field, modus_.impact_parameter(),
                      parameters_, projectile_target_interact_);
  // save evolution data
  if (!(modus_.is_box() && parameters_.outputclock->current_time() <
                               modus_.equilibration_time())) {
    for (const auto &output : outputs_) {
      if (output->is_dilepton_output() || output->is_photon_output() ||
          output->is_IC_output()) {
        continue;
      }

      output->at_intermediate_time(particles_, parameters_.outputclock,
                                   density_param_, event_info);

      // Thermodynamic output on the lattice versus time
      switch (dens_type_lattice_printout_) {
        case DensityType::Baryon:
          update_lattice(jmu_B_lat_.get(), lat_upd, DensityType::Baryon,
                         density_param_, particles_, false);
          output->thermodynamics_output(ThermodynamicQuantity::EckartDensity,
                                        DensityType::Baryon, *jmu_B_lat_);
          break;
        case DensityType::BaryonicIsospin:
          update_lattice(jmu_I3_lat_.get(), lat_upd,
                         DensityType::BaryonicIsospin, density_param_,
                         particles_, false);
          output->thermodynamics_output(ThermodynamicQuantity::EckartDensity,
                                        DensityType::BaryonicIsospin,
                                        *jmu_I3_lat_);
          break;
        case DensityType::None:
          break;
        default:
          update_lattice(jmu_custom_lat_.get(), lat_upd,
                         dens_type_lattice_printout_, density_param_,
                         particles_, false);
          output->thermodynamics_output(ThermodynamicQuantity::EckartDensity,
                                        dens_type_lattice_printout_,
                                        *jmu_custom_lat_);
      }
      if (printout_tmn_ || printout_tmn_landau_ || printout_v_landau_) {
        update_lattice(Tmn_.get(), lat_upd, dens_type_lattice_printout_,
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
}

template <typename Modus>
void Experiment<Modus>::update_potentials() {
  if (potentials_) {
    if (potentials_->use_symmetry() && jmu_I3_lat_ != nullptr) {
      update_lattice(jmu_I3_lat_.get(), LatticeUpdate::EveryTimestep,
                     DensityType::BaryonicIsospin, density_param_, particles_,
                     true);
    }
    if ((potentials_->use_skyrme() || potentials_->use_symmetry()) &&
        jmu_B_lat_ != nullptr) {
      update_lattice(jmu_B_lat_.get(), LatticeUpdate::EveryTimestep,
                     DensityType::Baryon, density_param_, particles_, true);
      const size_t UBlattice_size = UB_lat_->size();
      for (size_t i = 0; i < UBlattice_size; i++) {
        auto jB = (*jmu_B_lat_)[i];
        const FourVector flow_four_velocity_B =
            std::abs(jB.density()) > really_small ? jB.jmu_net() / jB.density()
                                                  : FourVector();
        double baryon_density = jB.density();
        ThreeVector baryon_grad_rho = jB.grad_rho();
        ThreeVector baryon_dj_dt = jB.dj_dt();
        ThreeVector baryon_rot_j = jB.rot_j();
        if (potentials_->use_skyrme()) {
          (*UB_lat_)[i] =
              flow_four_velocity_B * potentials_->skyrme_pot(baryon_density);
          (*FB_lat_)[i] = potentials_->skyrme_force(
              baryon_density, baryon_grad_rho, baryon_dj_dt, baryon_rot_j);
        }
        if (potentials_->use_symmetry() && jmu_I3_lat_ != nullptr) {
          auto jI3 = (*jmu_I3_lat_)[i];
          const FourVector flow_four_velocity_I3 =
              std::abs(jI3.density()) > really_small
                  ? jI3.jmu_net() / jI3.density()
                  : FourVector();
          (*UI3_lat_)[i] =
              flow_four_velocity_I3 *
              potentials_->symmetry_pot(jI3.density(), baryon_density);
          (*FI3_lat_)[i] = potentials_->symmetry_force(
              jI3.density(), jI3.grad_rho(), jI3.dj_dt(), jI3.rot_j(),
              baryon_density, baryon_grad_rho, baryon_dj_dt, baryon_rot_j);
        }
      }
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

    // Dileptons: shining of remaining resonances
    if (dilepton_finder_ != nullptr) {
      for (const auto &output : outputs_) {
        dilepton_finder_->shine_final(particles_, output.get(), true);
      }
    }
    // Find actions.
    for (const auto &finder : action_finders_) {
      actions.insert(finder->find_final_actions(particles_));
    }
    // Perform actions.
    while (!actions.is_empty()) {
      perform_action(*actions.pop(), particles_before_actions);
    }
    // loop until no more decays occur
  } while (interactions_total_ > interactions_old);

  // Dileptons: shining of stable particles at the end
  if (dilepton_finder_ != nullptr) {
    for (const auto &output : outputs_) {
      dilepton_finder_->shine_final(particles_, output.get(), false);
    }
  }
}

template <typename Modus>
void Experiment<Modus>::final_output(const int evt_num) {
  /* make sure the experiment actually ran (note: we should compare this
   * to the start time, but we don't know that. Therefore, we check that
   * the time is positive, which should heuristically be the same). */
  double E_mean_field = 0.0;
  if (likely(parameters_.labclock > 0)) {
    const uint64_t wall_actions_this_interval =
        wall_actions_total_ - previous_wall_actions_total_;
    const uint64_t interactions_this_interval = interactions_total_ -
                                                previous_interactions_total_ -
                                                wall_actions_this_interval;
    if (potentials_) {
      // using the lattice is necessary
      if ((jmu_B_lat_ != nullptr)) {
        E_mean_field =
            calculate_mean_field_energy(*potentials_, *jmu_B_lat_, parameters_);
      }
    }
    logg[LExperiment].info() << format_measurements(
        particles_, interactions_this_interval, conserved_initial_, time_start_,
        end_time_, E_mean_field, initial_mean_field_energy_);
    if (IC_output_switch_ && (particles_.size() == 0)) {
      // Verify there is no more energy in the system if all particles were
      // removed when crossing the hypersurface
      const double remaining_energy =
          conserved_initial_.momentum().x0() - total_energy_removed_;
      if (remaining_energy > really_small) {
        throw std::runtime_error(
            "There is remaining energy in the system although all particles "
            "were removed.\n"
            "E_remain = " +
            std::to_string(remaining_energy) + " [GeV]");
      } else {
        logg[LExperiment].info() << hline;
        logg[LExperiment].info()
            << "Time real: " << SystemClock::now() - time_start_;
        logg[LExperiment].info()
            << "Interactions before reaching hypersurface: "
            << interactions_total_ - wall_actions_total_ -
                   total_hypersurface_crossing_actions_;
        logg[LExperiment].info()
            << "Total number of particles removed on hypersurface: "
            << total_hypersurface_crossing_actions_;
      }
    } else {
      const double precent_discarded =
          interactions_total_ > 0
              ? static_cast<double>(discarded_interactions_total_) * 100.0 /
                    interactions_total_
              : 0.0;
      std::stringstream msg_discarded;
      msg_discarded
          << "Discarded interaction number: " << discarded_interactions_total_
          << " (" << precent_discarded
          << "% of the total interaction number including wall crossings)";

      logg[LExperiment].info() << hline;
      logg[LExperiment].info()
          << "Time real: " << SystemClock::now() - time_start_;
      logg[LExperiment].debug() << msg_discarded.str();

      if (parameters_.coll_crit == CollisionCriterion::Stochastic &&
          precent_discarded > 1.0) {
        // The choosen threshold of 1% is a heuristical value
        logg[LExperiment].warn()
            << msg_discarded.str()
            << "\nThe number of discarded interactions is large, which means "
               "the assumption for the stochastic criterion of\n1 interaction"
               "per particle per timestep is probably violated. Consider "
               "reducing the timestep size.";
      }

      logg[LExperiment].info() << "Final interaction number: "
                               << interactions_total_ - wall_actions_total_;
    }

    // Check if there are unformed particles
    int unformed_particles_count = 0;
    for (const auto &particle : particles_) {
      if (particle.formation_time() > end_time_) {
        unformed_particles_count++;
      }
    }
    if (unformed_particles_count > 0) {
      logg[LExperiment].warn(
          "End time might be too small. ", unformed_particles_count,
          " unformed particles were found at the end of the evolution.");
    }
  }

  auto event_info =
      fill_event_info(particles_, E_mean_field, modus_.impact_parameter(),
                      parameters_, projectile_target_interact_);

  for (const auto &output : outputs_) {
    output->at_eventend(particles_, evt_num, event_info);
  }
}

template <typename Modus>
void Experiment<Modus>::run() {
  const auto &mainlog = logg[LMain];
  for (int j = 0; j < nevents_; j++) {
    mainlog.info() << "Event " << j;

    // Sample initial particles, start clock, some printout and book-keeping
    initialize_new_event(j);
    /* In the ColliderModus, if the first collisions within the same nucleus are
     * forbidden, 'nucleon_has_interacted_', which records whether a nucleon has
     * collided with another nucleon, is initialized equal to false. If allowed,
     * 'nucleon_has_interacted' is initialized equal to true, which means these
     * incoming particles have experienced some fake scatterings, they can
     * therefore collide with each other later on since these collisions are not
     * "first" to them. */
    if (modus_.is_collider()) {
      if (!modus_.cll_in_nucleus()) {
        nucleon_has_interacted_.assign(modus_.total_N_number(), false);
      } else {
        nucleon_has_interacted_.assign(modus_.total_N_number(), true);
      }
    }
    /* In the ColliderModus, if Fermi motion is frozen, assign the beam momenta
     * to the nucleons in both the projectile and the target. */
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

    run_time_evolution();

    if (force_decays_) {
      do_final_decays();
    }

    // Output at event end
    final_output(j);
  }
}

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_EXPERIMENT_H_
