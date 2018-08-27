/*
 *    Copyright (c) 2013-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_EXPERIMENT_H_
#define SRC_INCLUDE_EXPERIMENT_H_

#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "actionfinderfactory.h"
#include "adaptiveparameters.h"
#include "chrono.h"
#include "decayactionsfinder.h"
#include "decayactionsfinderdilepton.h"
#include "energymomentumtensor.h"
#include "fourvector.h"
#include "grandcan_thermalizer.h"
#include "grid.h"
#include "outputparameters.h"
#include "pauliblocking.h"
#include "potentials.h"
#include "propagation.h"
#include "quantumnumbers.h"
#include "scatteractionphoton.h"
#include "scatteractionsfinder.h"
#include "thermalizationaction.h"
// Output
#include "binaryoutputcollisions.h"
#include "binaryoutputparticles.h"
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
  void initialize_new_event();

  /**
   * Runs the time evolution of an event with fixed-sized time steps,
   * adaptive time steps or without timesteps, from action to actions.
   * Within one timestep (fixed or adaptive) evolution from action to action
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
    return parameters_.outputclock.next_time();
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
  int n_fractional_photons_ = 100;

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

  /// Lattices for Skyme potentials
  std::unique_ptr<RectangularLattice<double>> UB_lat_ = nullptr;

  /// Lattices for symmetry potentials
  std::unique_ptr<RectangularLattice<double>> UI3_lat_ = nullptr;

  /// Lattices for the electric and magnetic components of the Skyme force
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

  /// system starting time of the simulation
  SystemTimePoint time_start_ = SystemClock::now();

  /// Type of density to be written to collision headers
  DensityType dens_type_ = DensityType::None;

  /// Pointer to additional parameters that are needed for adaptive time steps.
  std::unique_ptr<AdaptiveParameters> adaptive_parameters_ = nullptr;

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
void Experiment<Modus>::create_output(const std::string &format,
                                      const std::string &content,
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
    if (content == "Collisions" || content == "Dileptons" ||
        content == "Photons") {
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
    outputs_.emplace_back(make_unique<VtkOutput>(output_path, content));
  } else {
    log.error() << "Unknown combination of format (" << format
                << ") and content (" << content << "). Fix the config.";
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
 * The time after which the evolution is stopped. Note
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
 * \li \key None - No time steps are used. Cannot be used with potentials \n
 * \li \key Fixed - Fixed-sized time steps \n
 * \li \key Adaptive - Time steps with adaptive sizes.
 *
 * \key Metric_Type (string, optional, default = NoExpansion): \n
 * Select which kind of expansion the metric should have. This needs only be
 * specified for the sphere modus:
 * \li \key NoExpansion - Default SMASH run, with Minkowski metric \n
 * \li \key MasslessFRW - FRW expansion going as t^(1/2)
 * \li \key MassiveFRW - FRW expansion going as t^(2/3)
 * \li \key Exponential - FRW expansion going as e^(t/2)
 *
 * \key Expansion_Rate (double, optional, default = 0.1) \n
 * Corresponds to the speed of expansion of the universe in non minkowski
 * metrics if MetricType is any other than \key NoExpansion. \n
 * It corresponds to \f$b_r/l_0\f$ if the metric type is \key MasslessFRW or
 * \key MassiveFRW, and to the parameter b in the Exponential expansion where
 * \f$a(t) ~ e^{bt/2}\f$. \n
 *
 * \page input_collision_term_ Collision_Term
 *
 * \key Two_to_One (bool, optional, default = \key true) \n
 * Enable 2 <--> 1 processes (resonance formation and decays).
 *
 * \key Included_2to2 (list of 2 <--> 2 reactions,
 * optional, default = ["All"]) \n
 * List that contains all possible 2 <--> 2 process categories. Each process of
 * the listed category can be performed within the simulation. Possible
 * categories are: \li \key "Elastic" - elastic binary scatterings \li \key
 * "NN_to_NR" - nucleon + nucleon <--> nucleon + resonance \li \key "NN_to_DR" -
 * nucleon + nucleon <--> delta + resonance \li \key "KN_to_KN" - kaon + nucleon
 * <--> kaon + nucleon \li \key "KN_to_KDelta" - kaon + nucleon <--> kaon + dela
 * \li \key "Strangeness_exchange" - processes with strangeness exchange
 * \li \key "All" - include all binary processes, no necessity to list each
 * single category
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
 *  ρ → ππ and h₁(1170) → πρ, which gives 5 pions on average.
 * \li \key "strings" - Annihilation throug string fragmentation.
 *
 * \key Use_AQM (bool, optional, default = \key true) \n
 * Turn on AQM cross-sections for exotic combination of particles
 * (baryon-baryon cross-sections are scaled from proton-proton high energy
 * parametrization, for example). This includes both elastic and non-elastic
 * contributions; non-elastic contributions go through string fragmentation.
 * Turning off strings or elastic collisions while leaving this on will
 * result in the corresponding part of the AQM cross-sections to also be off.
 * Cross-sections parametrization are scaled according to
 * \f[ \frac{\sigma^{AQM}_{process}}{\sigma^{AQM}_{ref_process}}
 * \sigma^{param}_{ref_process}\f]
 * where \f$ \sigma^{AQM}_x = 40 \left( \frac{2}{3} \right)^{n_{meson}}
 * (1 - 0.4 x^s_1) (1 - 0.4 x^s_2) \f$, with \f$n_{meson}\f$ being the number
 * of mesons in the process, \f$x^s_{1,2}\f$ the fraction of strange quarks in
 * the participant. "process" is then a generic process and "ref_process" a
 * reference process such as PP for which solid parametrizations exist.
 * (\iref{Bass:1998ca})
 *
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
  bool no_coll = config.take({"Collision_Term", "No_Collisions"}, false);
  if ((parameters_.two_to_one || parameters_.included_2to2.any() ||
       parameters_.strings_switch) &&
      !no_coll) {
    auto scat_finder = make_unique<ScatterActionsFinder>(
        config, parameters_, nucleon_has_interacted_, modus_.total_N_number(),
        modus_.proj_N_number());
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
  ParticleData::formation_power_ =
      config.take({"Collision_Term", "Power_Particle_Formation"}, 1.);

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
   **/

  /*!\Userguide
   * \page input_general_
   *
   * \n
   * Example: Configuring General Properties
   * --------------
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
   * For the use of adaptive timesteps, change the \key Time_Step_Mode and
   * include the corresponding additional parameters:
   *\verbatim
       Time_Step_Mode: Adaptive
       Adaptive_Time_Step:
           Smoothing_Factor: 0.1
           Target_Missed_Actions: 0.01
           Allowed_Deviation: 2.5
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
   *                 computational frame or (optionally) only at the event end.
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
   *          and on a spatial lattice  versus time \ref Thermodynamics.
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
   *   - Used only for "Thermodynamics", see \ref Thermodynamics
   *
   * \note Output of coordinates for the "Collisions" content in
   *       the periodic box has a feature:
   *       \subpage collisions_output_in_box_modus_
   */

  /*!\Userguide
   * \page input_dileptons Dileptons
   * Enables dilepton output together with DecayActionsFinderDilepton.
   * By default, the extended OSCAR output is enabled. The dilepton output
   * format is identical to the collision output, it does however only contain
   * information about the dilepton decays at every timestep.
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
   * Enables the photon output together with the correpsoning
   * ScatterActionPhoton.
   *
   * Photons are treated perturbatively and are produced from binary
   * scattering processes. Their production follows the framework from Turbide
   * et al. described in \iref{Turbide:2006}. Following the perturbative
   * treatment, the produced photons do not contribute to the evolution of the
   * hadronic system. They are rather direcly printed to the photon output.
   * The mechanism for photon production is the following:
   * \li Look for hadronic interactions of particles that are also incoming
   * particles of a photon process. Currently, the latter include binary
   * scatterings of \f$ \pi \f$ and \f$ \rho \f$ mesons.
   * \li Perform the photon action and write the results to the photon output.
   * The final state particles are not of interested anymore as they are not
   * propagated further in the evolution. To account for the probabiity that
   * photon processes are significantly less likely than hadronic processes,
   * the produced photons are weighted according to the ratio of the photon
   * cross section to the hadronic cross section used the find the interaction.
   * This weight can be found in the weight element of the photon output.
   * \li Perform the original hadronic action based on which the photon action
   * was found. Propagate all final states particles throughout the hadronic
   * evolution as if no photon action had occured.
   *
   * As photons are produced very rarely, a lot of statistics is necessery to
   * yield useful results. Alternatively, it it possible to use fractional
   * photons. This means that for each produced photon, \f$ N_{\text{Frac}} \f$
   * photons are actually sampled with different kinematic properties so that
   * more phase space is covered. See \ref input_output_options_ on how to set
   *the flag.
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

  /* We can take away the Fermi motion flag, because the collider modus is
   * already initialized. We only need it when potentials are enabled, but we
   * always have to take it, otherwise SMASH will complain about unused
   * options. We have to provide a default value for modi other than Collider.
   */
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
   * For format of lattice output see \ref output_vtk_lattice_. To configure
   * output of the quantities on the lattice to vtk files see
   * \ref input_output_options_.
   *
   * \n
   * Examples: Configuring the Lattice
   * --------------
   * The following example configures the lattice with the origin in (0,0,0),
   * 20 cells of 10 fm size in each direction and with periodic boundary
   * conditions. The potential effects on the thresholds are taken into
   * consideration.
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
        UB_lat_ = make_unique<RectangularLattice<double>>(
            l, n, origin, periodic, LatticeUpdate::EveryTimestep);
        FB_lat_ = make_unique<
            RectangularLattice<std::pair<ThreeVector, ThreeVector>>>(
            l, n, origin, periodic, LatticeUpdate::EveryTimestep);
      }
      if (potentials_->use_symmetry()) {
        jmu_I3_lat_ = make_unique<DensityLattice>(l, n, origin, periodic,
                                                  LatticeUpdate::EveryTimestep);
        UI3_lat_ = make_unique<RectangularLattice<double>>(
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
    log.error(
        "If you want Thermodynamic VTK output, configure a lattice for it.");
  }

  // Store pointers to potential and lattice accessible for Action
  if (parameters_.potential_affect_threshold) {
    Action::input_potential(UB_lat_.get(), UI3_lat_.get(), potentials_.get());
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
const std::string hline(80, '-');

template <typename Modus>
void Experiment<Modus>::initialize_new_event() {
  const auto &log = logger<LogArea::Experiment>();

  random::set_seed(seed_);
  log.info() << "random number seed: " << seed_;
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

  // Reset the output clock
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
  // Print output headers
  log.info() << hline;
  log.info() << " Time       <Ediff>      <pdiff>  <scattrate>    <scatt>  "
                "<particles>   <timing>";
  log.info() << hline;
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
    /* Mark incoming nucleons as interacted - now they are permitted
     * to collide with nucleons from their native nucleus */
    for (const auto &incoming : action.incoming_particles()) {
      assert(incoming.id() >= 0);
      if (incoming.id() < modus_.total_N_number()) {
        nucleon_has_interacted_[incoming.id()] = true;
      }
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
  // Calculate Eckart rest frame density at the interaction point
  double rho = 0.0;
  if (dens_type_ != DensityType::None) {
    const FourVector r_interaction = action.get_interaction_point();
    constexpr bool compute_grad = false;
    rho = std::get<0>(rho_eckart(r_interaction.threevec(),
                                 particles_before_actions, density_param_,
                                 dens_type_, compute_grad));
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
  log.debug(~einhard::Green(), "✔ ", action);
  return true;
}

/**
 * Generate the tabulated string which will be printed to the screen when
 * SMASH is running
 *
 * \param[in] particles The interacting particles. Their information will be
 *            used to check the conservation of the total energy and momentum.
 *	      the total number of the particles will be used and printed as
 *	      well.
 * \param[in] scatterings_total Total number of the scatterings from the
 *            beginning to the current time step.
 * \param[in] scatterings_this_interval Number of the scatterings occur within
 *            the current timestep.
 * \param[in] conserved_initial Initial quantum numbers needed to check the
 *            conservations
 * \param[in] time_start Moment in the REAL WORLD when SMASH starts to run [s]
 * \param[in] time Current moment in SMASH [fm/c]
 * \return 'Current time in SMASH [fm/c]', 'Deviation of the energy from the
 *          initial value [GeV]', 'Deviation of the momentum from the initial
 *          value [GeV]', 'Averaged collisional rate [c/fm]', 'Number of the
 *          scatterings occur within the timestep', 'Total particle number',
 *          'Computing time consumed'
 */
std::string format_measurements(const Particles &particles,
                                uint64_t scatterings_total,
                                uint64_t scatterings_this_interval,
                                const QuantumNumbers &conserved_initial,
                                SystemTimePoint time_start, double time);

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

    // Perform forced thermalization if required
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
     *     compute new momenta according to equations of motion */
    if (potentials_) {
      update_potentials();
      update_momenta(&particles_, parameters_.labclock.timestep_duration(),
                     *potentials_, FB_lat_.get(), FI3_lat_.get());
    }

    /* (5) Expand universe if non-minkowskian metric; updates
     *     positions and momenta according to the selected expansion */
    if (metric_.mode_ != ExpansionMode::NoExpansion) {
      expand_space_time(&particles_, parameters_, metric_);
    }

    ++parameters_.labclock;

    /* (6) Check conservation laws.
     *
     * Check conservation of conserved quantities if potentials and string
     * fragmentation are off.  If potentials are on then momentum is conserved
     * only in average.  If string fragmentation is on, then energy and
     * momentum are only very roughly conserved in high-energy collisions. */
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
  /** \todo (Dima) is this part of the code still useful?
   *
   * if (thermalizer_) {
   *  thermalizer_->update_lattice(particles_, density_param_);
   *  thermalizer_->print_statistics(parameters_.labclock);
  }*/
  // save evolution data
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
                               density_param_, particles_, false);
        output->thermodynamics_output(ThermodynamicQuantity::EckartDensity,
                                      DensityType::Baryon, *jmu_B_lat_);
        break;
      case DensityType::BaryonicIsospin:
        update_density_lattice(jmu_I3_lat_.get(), lat_upd,
                               DensityType::BaryonicIsospin, density_param_,
                               particles_, false);
        output->thermodynamics_output(ThermodynamicQuantity::EckartDensity,
                                      DensityType::BaryonicIsospin,
                                      *jmu_I3_lat_);
        break;
      case DensityType::None:
        break;
      default:
        update_density_lattice(jmu_custom_lat_.get(), lat_upd,
                               dens_type_lattice_printout_, density_param_,
                               particles_, false);
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
                             DensityType::Baryon, density_param_, particles_,
                             true);
      const size_t UBlattice_size = UB_lat_->size();
      for (size_t i = 0; i < UBlattice_size; i++) {
        (*UB_lat_)[i] = potentials_->skyrme_pot((*jmu_B_lat_)[i].density());
        (*FB_lat_)[i] = potentials_->skyrme_force(
            (*jmu_B_lat_)[i].density(), (*jmu_B_lat_)[i].grad_rho(),
            (*jmu_B_lat_)[i].dj_dt(), (*jmu_B_lat_)[i].rot_j());
      }
    }
    if (potentials_->use_symmetry() && jmu_I3_lat_ != nullptr) {
      update_density_lattice(jmu_I3_lat_.get(), LatticeUpdate::EveryTimestep,
                             DensityType::BaryonicIsospin, density_param_,
                             particles_, true);
      const size_t UI3lattice_size = UI3_lat_->size();
      for (size_t i = 0; i < UI3lattice_size; i++) {
        (*UI3_lat_)[i] = potentials_->symmetry_pot((*jmu_I3_lat_)[i].density());
        (*FI3_lat_)[i] = potentials_->symmetry_force(
            (*jmu_I3_lat_)[i].grad_rho(), (*jmu_I3_lat_)[i].dj_dt(),
            (*jmu_I3_lat_)[i].rot_j());
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
  const auto &log = logger<LogArea::Experiment>();
  /* make sure the experiment actually ran (note: we should compare this
   * to the start time, but we don't know that. Therefore, we check that
   * the time is positive, which should heuristically be the same). */
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
    // if there are no particles no interactions happened
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
    for (const auto &particle : particles_) {
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

    // Sample initial particles, start clock, some printout and book-keeping
    initialize_new_event();
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

    // Output at event start
    for (const auto &output : outputs_) {
      output->at_eventstart(particles_, j);
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

#endif  // SRC_INCLUDE_EXPERIMENT_H_
