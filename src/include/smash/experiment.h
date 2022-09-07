/*
 *    Copyright (c) 2013-2022
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
#include "fields.h"
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
#ifdef SMASH_USE_RIVET
#include "rivetoutput.h"
#endif
#include "icoutput.h"
#include "oscaroutput.h"
#include "thermodynamiclatticeoutput.h"
#include "thermodynamicoutput.h"
#ifdef SMASH_USE_ROOT
#include "rootoutput.h"
#endif
#include "freeforallaction.h"
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

  /**
   * \ingroup exception
   * Exception class that is thrown if the requested output path in the
   * Experiment factory is not existing.
   */
  struct NonExistingOutputPathRequest : public std::invalid_argument {
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
   * Runs the time evolution of an event with fixed-sized time steps or without
   * timesteps, from action to actions.
   * Within one timestep (fixed) evolution from action to action
   * is invoked.
   *
   * \param[in] t_end time until run_time_evolution is run, in SMASH this is the
   *                  configured end_time, but it might differ if SMASH is used
   *                  as an external library
   */
  void run_time_evolution(const double t_end, ActionList extra_actions);

  /**
   * Performs the final decays of an event
   *
   * \throws runtime_error if found actions cannot be performed
   */
  void do_final_decays();

  /// Output at the end of an event
  void final_output();

  /**
   * Provides external access to SMASH particles. This is helpful if SMASH
   * is used as a 3rd-party library.
   */
  // todo(oliiny): this should be made compatible with JetScape on
  // the JetScape side
  Particles *first_ensemble() { return &ensembles_[0]; }
  /// Getter for all ensembles
  std::vector<Particles> *all_ensembles() { return &ensembles_; }

  /**
   * Provides external access to SMASH calculation modus. This is helpful if
   * SMASH is used as a 3rd-party library.
   */
  Modus *modus() { return &modus_; }

 private:
  /**
   * Perform the given action.
   *
   * \param[in] action The action to perform
   * \param[in] i_ensemble index of ensemble in which action is performed
   * \param[in] include_pauli_blocking wheter to take Pauli blocking into
   *                                   account. Skipping Pauli blocking is
   *                                   useful for example for final decays.
   * \return False if the action is
   *                 rejected either due to invalidity or
   *                 Pauli-blocking, or true if it's accepted and performed.
   */
  bool perform_action(Action &action, int i_ensemble,
                      bool include_pauli_blocking = true);
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
   * \param[in, out] particles Particles to be propagated
   */
  void propagate_and_shine(double to_time, Particles &particles);

  /**
   * Performs all the propagations and actions during a certain time interval
   * neglecting the influence of the potentials. This function is called in
   * either the time stepless cases or the cases with time steps.
   *
   * \param[in, out] actions Actions occur during a certain time interval.
   *                 They provide the ending times of the propagations and
   *                 are updated during the time interval.
   * \param[in]      i_ensemble index of ensemble to be evolved
   * \param[in]      end_time_propagation time until propagation should be
   *                 performed
   * \param[in]      end_time_run time until the whole evolution is run
   */
  void run_time_evolution_timestepless(Actions &actions, int i_ensemble,
                                       const double end_time_propagation,
                                       const double end_time_run);

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
    if (parameters_.coll_crit == CollisionCriterion::Stochastic) {
      return parameters_.fixed_min_cell_length;
    }
    return std::sqrt(4 * dt * dt + max_transverse_distance_sqr_);
  }

  /// Shortcut for next output time
  double next_output_time() const {
    return parameters_.outputclock->next_time();
  }

  /**
   * Counts the number of ensembles in wich interactions took place at the end
   * of an event
   */
  void count_nonempty_ensembles();

  /**
   * Checks wether the desired number events have been calculated
   *
   * \return wether the experiment is is_finished
   */
  bool is_finished();

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

  /// Complete particle list, all ensembles in one vector
  std::vector<Particles> ensembles_;

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
   * Whether the projectile and the target collided.
   * One value for each ensemble.
   */
  std::vector<bool> projectile_target_interact_;

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

  /// 4-current for j_QBS lattice output
  std::unique_ptr<DensityLattice> j_QBS_lat_;

  /// Baryon density on the lattice
  std::unique_ptr<DensityLattice> jmu_B_lat_;

  /// Isospin projection density on the lattice
  std::unique_ptr<DensityLattice> jmu_I3_lat_;

  /// Electric charge density on the lattice
  std::unique_ptr<DensityLattice> jmu_el_lat_;

  /// Mean-field A^mu on the lattice
  std::unique_ptr<FieldsLattice> fields_lat_;

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
   * Lattices for Skyrme or VDF potentials (evaluated in the local rest frame)
   * times the baryon flow 4-velocity
   */
  std::unique_ptr<RectangularLattice<FourVector>> UB_lat_ = nullptr;

  /**
   * Lattices for symmetry potentials (evaluated in the local rest frame) times
   * the isospin flow 4-velocity
   */
  std::unique_ptr<RectangularLattice<FourVector>> UI3_lat_ = nullptr;

  /**
   * Lattices for the electric and magnetic components of the Skyrme or VDF
   * force
   */
  std::unique_ptr<RectangularLattice<std::pair<ThreeVector, ThreeVector>>>
      FB_lat_;

  /// Lattices for the electric and magnetic component of the symmetry force
  std::unique_ptr<RectangularLattice<std::pair<ThreeVector, ThreeVector>>>
      FI3_lat_;

  /// Lattices for electric and magnetic field in fm^-2
  std::unique_ptr<RectangularLattice<std::pair<ThreeVector, ThreeVector>>>
      EM_lat_;

  /// Lattices of energy-momentum tensors for printout
  std::unique_ptr<RectangularLattice<EnergyMomentumTensor>> Tmn_;

  /// Auxiliary lattice for values of jmu at a time step t0
  std::unique_ptr<RectangularLattice<FourVector>> old_jmu_auxiliary_;
  /// Auxiliary lattice for values of jmu at a time step t0 + dt
  std::unique_ptr<RectangularLattice<FourVector>> new_jmu_auxiliary_;
  /// Auxiliary lattice for calculating the four-gradient of jmu
  std::unique_ptr<RectangularLattice<std::array<FourVector, 4>>>
      four_gradient_auxiliary_;

  /// Auxiliary lattice for values of Amu at a time step t0
  std::unique_ptr<RectangularLattice<FourVector>> old_fields_auxiliary_;
  /// Auxiliary lattice for values of Amu at a time step t0 + dt
  std::unique_ptr<RectangularLattice<FourVector>> new_fields_auxiliary_;
  /// Auxiliary lattice for calculating the four-gradient of Amu
  std::unique_ptr<RectangularLattice<std::array<FourVector, 4>>>
      fields_four_gradient_auxiliary_;

  /// Whether to print the energy-momentum tensor
  bool printout_tmn_ = false;

  /// Whether to print the energy-momentum tensor in Landau frame
  bool printout_tmn_landau_ = false;

  /// Whether to print the 4-velocity in Landau frame
  bool printout_v_landau_ = false;

  /// Whether to print the Q, B, S 4-currents
  bool printout_j_QBS_ = false;

  /// Whether to print the thermodynamics quantities evaluated on the lattices
  bool printout_lattice_td_ = false;

  /// Whether to print the thermodynamics quantities evaluated on the lattices,
  /// point by point, in ASCII format
  bool printout_full_lattice_ascii_td_ = false;

  /// Whether to print the thermodynamics quantities evaluated on the lattices,
  /// point by point, in Binary format
  bool printout_full_lattice_binary_td_ = false;

  /// Whether to print the thermodynamics quantities evaluated on the lattices,
  /// point by point, in any format
  bool printout_full_lattice_any_td_ = false;

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
  const int nevents_ = 0;

  /**
   * The number of ensembles, in which interactions take place, to be
   * calculated.
   *
   * Can be specified as an inout instead of the number of events. In
   * this case events will be calculated until this number of ensembles
   * is reached.
   */
  int minimum_nonempty_ensembles_ = 0;

  /**
   * The way in which the number of calculated events is specified.
   *
   * Can be either a fixed number of simulated events or a minimum number
   * of events that contain interactions.
   */
  EventCounting event_counting_ = EventCounting::Invalid;

  /// Current event
  int event_ = 0;

  /// Number of ensembles containing an interaction
  int nonempty_ensembles_ = 0;

  /**
   * Maximum number of events to be calculated in order obtain the desired
   * number of non-empty events using the MinimumNonemptyEnsembles option.
   */
  int max_events_ = 0;

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

  /**
   * Maximal distance at which particles can interact in case of the geometric
   * criterion, squared
   */
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

  /// This indicates whether kinematic cuts are enabled for the IC output
  bool kinematic_cuts_for_IC_output_ = false;

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
  } else if (content == "Thermodynamics" &&
             (format == "Lattice_ASCII" || format == "Lattice_Binary")) {
    printout_full_lattice_any_td_ = true;
    if (format == "Lattice_ASCII") {
      printout_full_lattice_ascii_td_ = true;
    }
    if (format == "Lattice_Binary") {
      printout_full_lattice_binary_td_ = true;
    }
    outputs_.emplace_back(make_unique<ThermodynamicLatticeOutput>(
        output_path, content, out_par, printout_full_lattice_ascii_td_,
        printout_full_lattice_binary_td_));
  } else if (content == "Thermodynamics" && format == "VTK") {
    printout_lattice_td_ = true;
    outputs_.emplace_back(
        make_unique<VtkOutput>(output_path, content, out_par));
  } else if (content == "Initial_Conditions" && format == "ASCII") {
    outputs_.emplace_back(
        make_unique<ICOutput>(output_path, "SMASH_IC", out_par));
  } else if ((format == "HepMC") || (format == "HepMC_asciiv3") ||
             (format == "HepMC_treeroot")) {
#ifdef SMASH_USE_HEPMC
    if (content == "Particles") {
      if ((format == "HepMC") || (format == "HepMC_asciiv3")) {
        outputs_.emplace_back(make_unique<HepMcOutput>(
            output_path, "SMASH_HepMC_particles", false, "asciiv3"));
      } else if (format == "HepMC_treeroot") {
#ifdef SMASH_USE_HEPMC_ROOTIO
        outputs_.emplace_back(make_unique<HepMcOutput>(
            output_path, "SMASH_HepMC_particles", false, "root"));
#else
        logg[LExperiment].error(
            "Requested HepMC_treeroot output not available, "
            "ROOT or HepMC3 ROOTIO missing or not found by cmake.");
#endif
      }
    } else if (content == "Collisions") {
      if ((format == "HepMC") || (format == "HepMC_asciiv3")) {
        outputs_.emplace_back(make_unique<HepMcOutput>(
            output_path, "SMASH_HepMC_collisions", true, "asciiv3"));
      } else if (format == "HepMC_treeroot") {
#ifdef SMASH_USE_HEPMC_ROOTIO
        outputs_.emplace_back(make_unique<HepMcOutput>(
            output_path, "SMASH_HepMC_collisions", true, "root"));
#else
        logg[LExperiment].error(
            "Requested HepMC_treeroot output not available, "
            "ROOT or HepMC3 ROOTIO missing or not found by cmake.");
#endif
      }
    } else {
      logg[LExperiment].error(
          "HepMC only available for Particles and "
          "Collisions content. Requested for " +
          content + ".");
    }
#else
    logg[LExperiment].error(
        "HepMC output requested, but HepMC support not compiled in");
#endif
  } else if (content == "Coulomb" && format == "VTK") {
    outputs_.emplace_back(
        make_unique<VtkOutput>(output_path, "Fields", out_par));
  } else if (content == "Rivet") {
#ifdef SMASH_USE_RIVET
    // flag to ensure that the Rivet format has not been already assigned
    static bool rivet_format_already_selected = false;
    // if the next check is true, then we are trying to assign the format twice
    if (rivet_format_already_selected) {
      logg[LExperiment].warn(
          "Rivet output format can only be one, either YODA or YODA-full. "
          "Only your first valid choice will be used.");
      return;
    }
    if (format == "YODA") {
      outputs_.emplace_back(
          make_unique<RivetOutput>(output_path, "SMASH_Rivet", false, out_par));
      rivet_format_already_selected = true;
    } else if (format == "YODA-full") {
      outputs_.emplace_back(make_unique<RivetOutput>(
          output_path, "SMASH_Rivet_full", true, out_par));
      rivet_format_already_selected = true;
    } else {
      logg[LExperiment].error("Rivet format " + format +
                              "not one of YODA or YODA-full");
    }
#else
    logg[LExperiment].error(
        "Rivet output requested, but Rivet support not compiled in");
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
 * \key Nevents (int, optional, no default): \n
 * Number of events to calculate.\n
 * One can instead define a minimum number of ensembles that contain
 * interactions, see \subpage minimum_nonempty_ensembles.
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
 * in particles.txt as well. Do not uncomment d' or make sure to exclude 2-body
 * reactions involving the d' (i.e. no \key "PiDeuteron_to_pidprime" and
 * \key "NDeuteron_to_Ndprime" in \key Included_2to2). Otherwise, the deuteron
 * reactions are implicitly double-counted.
 * \li \key "A3_Nuclei_4to2" - Create or destroy
 * A = 3 nuclei (triton, He-3, hypertriton) by \f$4 \leftrightarrow 2\f$
 * catalysis reactions such as \f$ X NNN \leftrightarrow X t \f$, where
 * X can be a pion, nucleon, or antinucleon.
 * \li \key "NNbar_5to2" - 5-to-2 backreaction for NNbar  annihilation:
 * \f$ \pi^0\pi^+\pi^-\pi^+\pi^- \rightarrow N\bar{N}\f$. Since detailed balance
 * is enforced, NNbar_Treatment has to be set to "two to five" for this option.
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
 * \li \key "two to five" - Direct Annhilation of NNbar to 5π, matching the
 *  resonance treatment: \f$ N\bar{N}\rightarrow\pi^0\pi^+\pi^-\pi^+\pi^-\f$ .
 *  This option requires "NNbar_5to2" to be enabled in Multi_Particle_Reactions.
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
 * \page minimum_nonempty_ensembles Minimum non-empty Ensembles
 * Instead of defining the number of Events it is possible to define a minimum
 * number
 * of ensembles in which an interaction took place. Using this option, events
 * will
 * be calculated until the desired number of non-empty ensembles is generated.
 * \n
 * number
 * of ensembles is equal to the number of events, so that this option will
 * provide
 * the desired number of non-empty events.\n
 *
 * \key Number (int, required): \n
 * The number of desired non-empty ensembles\n
 *
 * \key Maximum_Ensembles_Run (int, requiured): \n
 * Maximum number of ensembles run. This number serves as a safeguard
 * against SMASH unexpectedly running for a long time.
 *
 * **Examples**\n
 * In the following example, the number of desired non-empty events is 1000
 * with a maximum number of 2000 events to be calculated. In this case the
 * calculation will stop either if 1000 events are not empty or 2000 events
 * have been calculated.
 * \verbatim
  General:
   Modus: Collider
    Minimum_Nonempty_Ensembles:
      Number: 1000
      Maximum_Ensembles_Run: 2000
    Ensembles: 1
  \endverbatim
 * In contrast to the first example, we use 20 parallel ensembles here.
 * Here, the maximum number of ensembles run is 2000.
 * The calculation will continue until either this number of ensembles
 * is reached or 1000 ensembles contain interactions. Note that an event
 * consists of 20 ensembles. The 20 ensembles run in parallel, so the number of
 * non-empty ensembles in the ouput is between 1000 and 1019.
 * \verbatim
  General:
    Modus: Collider
    Minimum_Nonempty_Ensembles:
      Number: 1000
      Maximum_Ensembles_Run: 2000
    Ensembles: 20
  \endverbatim
 */
template <typename Modus>
Experiment<Modus>::Experiment(Configuration config, const bf::path &output_path)
    : parameters_(create_experiment_parameters(config)),
      density_param_(DensityParameters(parameters_)),
      modus_(config["Modi"], parameters_),
      ensembles_(parameters_.n_ensembles),
      nevents_(config.take({"General", "Nevents"}, 0)),
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

  if (config.has_value({"General", "Minimum_Nonempty_Ensembles"})) {
    if (config.has_value({"General", "Nevents"})) {
      throw std::invalid_argument(
          "Please specify either Nevents or Minimum_Nonempty_Ensembles.");
    }
    event_counting_ = EventCounting::MinimumNonEmpty;
    minimum_nonempty_ensembles_ =
        config.take({"General", "Minimum_Nonempty_Ensembles", "Number"});
    int max_ensembles = config.take(
        {"General", "Minimum_Nonempty_Ensembles", "Maximum_Ensembles_Run"});
    max_events_ =
        std::ceil(static_cast<double>(max_ensembles) / parameters_.n_ensembles);
  } else {
    event_counting_ = EventCounting::FixedNumber;
  }

  // covariant derivatives can only be done with covariant smearing
  if (parameters_.derivatives_mode == DerivativesMode::CovariantGaussian &&
      parameters_.smearing_mode != SmearingMode::CovariantGaussian) {
    throw std::invalid_argument(
        "Covariant Gaussian derivatives only make sense for Covariant Gaussian "
        "smearing!");
  }

  // for triangular smearing:
  // the weight needs to be larger than 1./7. for the center cell to contribute
  // more than the surrounding cells
  if (parameters_.smearing_mode == SmearingMode::Discrete &&
      parameters_.discrete_weight < (1. / 7.)) {
    throw std::invalid_argument(
        "The central weight for discrete smearing should be >= 1./7.");
  }

  if (parameters_.coll_crit == CollisionCriterion::Stochastic &&
      (time_step_mode_ != TimeStepMode::Fixed || !use_grid_)) {
    throw std::invalid_argument(
        "The stochastic criterion can only be employed for fixed time step "
        "mode and with a grid!");
  }

  logg[LExperiment].info("Using ", parameters_.testparticles,
                         " testparticles per particle.");
  logg[LExperiment].info("Using ", parameters_.n_ensembles,
                         " parallel ensembles.");

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
    action_finders_.emplace_back(make_unique<DecayActionsFinder>(
        parameters_.res_lifetime_factor, parameters_.do_weak_decays));
  }
  bool no_coll = config.take({"Collision_Term", "No_Collisions"}, false);
  if ((parameters_.two_to_one || parameters_.included_2to2.any() ||
       parameters_.included_multi.any() || parameters_.strings_switch) &&
      !no_coll) {
    auto scat_finder = make_unique<ScatterActionsFinder>(config, parameters_);
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

    double rapidity_cut = 0.0;
    if (config.has_value({"Output", "Initial_Conditions", "Rapidity_Cut"})) {
      rapidity_cut =
          config.take({"Output", "Initial_Conditions", "Rapidity_Cut"});
      if (rapidity_cut <= 0.0) {
        logg[LInitialConditions].fatal()
            << "Rapidity cut for initial conditions configured as abs(y) < "
            << rapidity_cut << " is unreasonable. \nPlease choose a positive, "
            << "non-zero value or employ SMASH without rapidity cut.";
        throw std::runtime_error(
            "Kinematic cut for initial conditions malconfigured.");
      }
    }

    if (modus_.calculation_frame_is_fixed_target() && rapidity_cut != 0.0) {
      throw std::runtime_error(
          "Rapidity cut for initial conditions output is not implemented "
          "in the fixed target calculation frame. \nPlease use "
          "\"center of velocity\" or \"center of mass\" as a "
          "\"Calculation_Frame\" instead.");
    }

    double transverse_momentum_cut = 0.0;
    if (config.has_value({"Output", "Initial_Conditions", "pT_Cut"})) {
      transverse_momentum_cut =
          config.take({"Output", "Initial_Conditions", "pT_Cut"});
      if (transverse_momentum_cut <= 0.0) {
        logg[LInitialConditions].fatal()
            << "transverse momentum cut for initial conditions configured as "
               "pT < "
            << rapidity_cut << " is unreasonable. \nPlease choose a positive, "
            << "non-zero value or employ SMASH without pT cut.";
        throw std::runtime_error(
            "Kinematic cut for initial conditions malconfigured.");
      }
    }

    if (rapidity_cut > 0.0 || transverse_momentum_cut > 0.0) {
      kinematic_cuts_for_IC_output_ = true;
    }

    if (rapidity_cut > 0.0 && transverse_momentum_cut > 0.0) {
      logg[LInitialConditions].info()
          << "Extracting initial conditions in kinematic range: "
          << -rapidity_cut << " <= y <= " << rapidity_cut
          << "; pT <= " << transverse_momentum_cut << " GeV.";
    } else if (rapidity_cut > 0.0) {
      logg[LInitialConditions].info()
          << "Extracting initial conditions in kinematic range: "
          << -rapidity_cut << " <= y <= " << rapidity_cut << ".";
    } else if (transverse_momentum_cut > 0.0) {
      logg[LInitialConditions].info()
          << "Extracting initial conditions in kinematic range: pT <= "
          << transverse_momentum_cut << " GeV.";
    } else {
      logg[LInitialConditions].info()
          << "Extracting initial conditions without kinematic cuts.";
    }

    action_finders_.emplace_back(make_unique<HyperSurfaceCrossActionsFinder>(
        proper_time, rapidity_cut, transverse_momentum_cut));
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
  logg[LExperiment].trace(SMASH_SOURCE_LOCATION,
                          " create OutputInterface objects");

  auto output_conf = config["Output"];
  /*!\Userguide
   * \page output_general_ Output
   *
   * \section output_directory_ Output directory
   *
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
   * \section output_contents_ Output content
   *
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
   *      \ref format_binary_, \ref format_root, \ref format_vtk, \ref
   * output_hepmc_
   * - \b Collisions List of interactions: collisions, decays, box wall
   *                 crossings and forced thermalizations. Information about
   *                 incoming, outgoing particles and the interaction itself
   *                 is printed out.
   *   - Available formats: \ref format_oscar_collisions, \ref format_binary_,
   *                 \ref format_root, \subpage output_hepmc_
   * - \b Dileptons  Special dilepton output, see \subpage output_dileptons.
   *   - Available formats: \ref format_oscar_collisions,
   *                   \ref format_binary_ and \ref format_root
   * - \b Photons    Special photon output, see \subpage output_photons.
   *   - Available formats: \ref format_oscar_collisions,
   *                   \ref format_binary_ and \ref format_root.
   * - \b Thermodynamics   This output allows to print out thermodynamic
   *          quantities, see \ref Thermodynamics.
   *    - Available formats: \ref thermodyn_output_user_guide_,
   *      \ref thermodyn_lattice_output_,
   *      \ref output_vtk_lattice_
   * - \b Initial_Conditions  Special initial conditions output, see
   *                          \subpage input_ic for details
   *   - Available formats: \ref format_oscar_particlelist, \ref
   * IC_output_user_guide_
   * - \b Rivet Run Rivet analysis on generated events and output
   *    results, see \subpage rivet_output_user_guide_ for details.
   *    - Available formats: \ref rivet_output_user_guide_
   *
   *
   * \n
   *
   * \section list_of_output_formats Output formats
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
   *   - Used for "Thermodynamics" and "Initial_Conditions", see
   * \subpage thermodyn_output_user_guide_
   * \subpage thermodyn_lattice_output_
   * \subpage IC_output_user_guide_
   * - \b "HepMC_asciiv3", \b "HepMC_treeroot" - HepMC3 human-readble asciiv3 or
   *   Tree ROOT format see \ref output_hepmc_ for details
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
    // potentials need density calculation parameters from parameters_
    potentials_ = make_unique<Potentials>(config["Potentials"], parameters_);
    // make sure that vdf potentials are not used together with Skyrme
    // or symmetry potentials
    if (potentials_->use_skyrme() && potentials_->use_vdf()) {
      throw std::runtime_error(
          "Can't use Skyrme and VDF potentials at the same time!");
    }
    if (potentials_->use_symmetry() && potentials_->use_vdf()) {
      throw std::runtime_error(
          "Can't use symmetry and VDF potentials at the same time!");
    }
    if (potentials_->use_skyrme()) {
      logg[LExperiment].info() << "Skyrme potentials are:\n";
      logg[LExperiment].info()
          << "\t\tSkyrme_A [MeV] = " << potentials_->skyrme_a() << "\n";
      logg[LExperiment].info()
          << "\t\tSkyrme_B [MeV] = " << potentials_->skyrme_b() << "\n";
      logg[LExperiment].info()
          << "\t\t    Skyrme_tau = " << potentials_->skyrme_tau() << "\n";
    }
    if (potentials_->use_symmetry()) {
      logg[LExperiment].info()
          << "Symmetry potential is:"
          << "\n   S_pot [MeV] = " << potentials_->symmetry_S_pot() << "\n";
    }
    if (potentials_->use_vdf()) {
      logg[LExperiment].info() << "VDF potential parameters are:\n";
      logg[LExperiment].info() << "\t\tsaturation density [fm^-3] = "
                               << potentials_->saturation_density() << "\n";
      for (int i = 0; i < potentials_->number_of_terms(); i++) {
        logg[LExperiment].info()
            << "\t\tCoefficient_" << i + 1 << " = "
            << 1000.0 * (potentials_->coeffs())[i] << " [MeV]   \t Power_"
            << i + 1 << " = " << (potentials_->powers())[i] << "\n";
      }
    }
    // if potentials are on, derivatives need to be calculated
    if (parameters_.derivatives_mode == DerivativesMode::Off &&
        parameters_.field_derivatives_mode == FieldDerivativesMode::ChainRule) {
      throw std::invalid_argument(
          "Derivatives are necessary for running with potentials.\n"
          "Derivatives_Mode: \"Off\" only makes sense for "
          "Field_Derivatives_Mode: \"Direct\"!\nUse \"Covariant Gaussian\" or "
          "\"Finite difference\".");
    }
    // for computational efficiency, we want to turn off the derivatives of jmu
    // and the rest frame density derivatives if direct derivatives are used
    if (parameters_.field_derivatives_mode == FieldDerivativesMode::Direct) {
      parameters_.derivatives_mode = DerivativesMode::Off;
      parameters_.rho_derivatives_mode = RestFrameDensityDerivativesMode::Off;
    }
    switch (parameters_.derivatives_mode) {
      case DerivativesMode::CovariantGaussian:
        logg[LExperiment].info() << "Covariant Gaussian derivatives are ON";
        break;
      case DerivativesMode::FiniteDifference:
        logg[LExperiment].info() << "Finite difference derivatives are ON";
        break;
      case DerivativesMode::Off:
        logg[LExperiment].info() << "Gradients of baryon current are OFF";
        break;
    }
    switch (parameters_.rho_derivatives_mode) {
      case RestFrameDensityDerivativesMode::On:
        logg[LExperiment].info() << "Rest frame density derivatives are ON";
        break;
      case RestFrameDensityDerivativesMode::Off:
        logg[LExperiment].info() << "Rest frame density derivatives are OFF";
        break;
    }
    // direct or chain rule derivatives only make sense for the VDF potentials
    if (potentials_->use_vdf()) {
      switch (parameters_.field_derivatives_mode) {
        case FieldDerivativesMode::ChainRule:
          logg[LExperiment].info() << "Chain rule field derivatives are ON";
          break;
        case FieldDerivativesMode::Direct:
          logg[LExperiment].info() << "Direct field derivatives are ON";
          break;
      }
    }
    /*
     * Necessary safety checks
     */
    // VDF potentials need derivatives of rest frame density or fields
    if (potentials_->use_vdf() && (parameters_.rho_derivatives_mode ==
                                       RestFrameDensityDerivativesMode::Off &&
                                   parameters_.field_derivatives_mode ==
                                       FieldDerivativesMode::ChainRule)) {
      throw std::runtime_error(
          "Can't use VDF potentials without rest frame density derivatives or "
          "direct field derivatives!");
    }
    // potentials require using gradients
    if (parameters_.derivatives_mode == DerivativesMode::Off &&
        parameters_.field_derivatives_mode == FieldDerivativesMode::ChainRule) {
      throw std::runtime_error(
          "Can't use potentials without gradients of baryon current (Skyrme, "
          "VDF)"
          " or direct field derivatives (VDF)!");
    }
    // direct field derivatives only make sense for the VDF potentials
    if (!(potentials_->use_vdf()) &&
        parameters_.field_derivatives_mode == FieldDerivativesMode::Direct) {
      throw std::invalid_argument(
          "Field_Derivatives_Mode: \"Direct\" only makes sense for the VDF "
          "potentials!\nUse Field_Derivatives_Mode: \"Chain Rule\" or comment "
          "this option out (Chain Rule is default)");
    }
  }

  // information about the type of smearing
  switch (parameters_.smearing_mode) {
    case SmearingMode::CovariantGaussian:
      logg[LExperiment].info() << "Smearing type: Covariant Gaussian";
      break;
    case SmearingMode::Discrete:
      logg[LExperiment].info() << "Smearing type: Discrete with weight = "
                               << parameters_.discrete_weight;
      break;
    case SmearingMode::Triangular:
      logg[LExperiment].info() << "Smearing type: Triangular with range = "
                               << parameters_.triangular_range;
      break;
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
   * \ref output_vtk_lattice_), the Thermodynamic Lattice Output (see
   * \ref thermodyn_lattice_output_ ) or the \key Potentials_Affect_Thresholds
   * option
   * is enabled. \n
   * The following parameters are only required, if the \key Lattice section is
   * used in the configuration file. Otherwise, no lattice will be used at all.
   *
   *
   * \key Sizes (array<double,3>, optional, default depends on modus): \n
   *      Sizes of lattice in x, y, z directions in fm.
   *
   * \key Cell_Number (array<int,3>, optional, default depends on modus): \n
   *      Number of cells in x, y, z directions.
   *
   * \key Origin (array<double,3>, optional, default depends on modus): \n
   *      Coordinates of the left, down, near corner of the lattice in fm.
   *
   * \key Periodic (bool, optional, default true for Box modus,
   *                                        false for other modi): \n
   *      Use periodic continuation or not. With periodic continuation
   *      x + i * lx is equivalent to x, same for y, z.
   *
   * \key Potentials_Affect_Thresholds (bool, optional, default = false): \n
   * Include potential effects, since mean field potentials change the threshold
   * energies of the actions.
   *
   * For information on the format of the lattice output see
   * \ref output_vtk_lattice_ or \ref thermodyn_lattice_output_. To configure
   * the thermodynamic output, see \ref input_output_options_.
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
   *\n
   * In case of Collider, Box, and Sphere modus
   * (see input_general_ for choosing modus)
   * there is also an option to set up lattice automatically.
   * For example, for Collider modus
   *
   *\verbatim
   Lattice:
   \endverbatim
   * sets up a lattice that (heuristically) covers causally allowed particle
   * positions until the end time of the simulation
   * (see input_general_ for choosing end time).
   * The lattice may also be automatically contracted in z direction depending
   * on the chosen way of density calculation.
   *
   * Another example for Box modus:
   *\verbatim
   Lattice:
       Cell_Number:    [20, 20, 20]
   \endverbatim
   * sets up a periodic lattice that matches box sizes.
   */

  // Create lattices
  if (config.has_value({"Lattice"})) {
    std::array<double, 3> l_default{20., 20., 20.};
    std::array<int, 3> n_default{10, 10, 10};
    std::array<double, 3> origin_default{-20., -20., -20.};
    bool periodic_default = false;
    if (modus_.is_collider()) {
      // Estimates on how far particles could get in x, y, z
      const double gam = modus_.sqrt_s_NN() / (2.0 * nucleon_mass);
      const double max_z = 5.0 / gam + end_time_;
      const double estimated_max_transverse_velocity = 0.7;
      const double max_xy = 5.0 + estimated_max_transverse_velocity * end_time_;
      origin_default = {-max_xy, -max_xy, -max_z};
      l_default = {2 * max_xy, 2 * max_xy, 2 * max_z};
      // Go for approximately 0.8 fm cell size and contract
      // lattice in z by gamma factor
      const int n_xy = std::ceil(2 * max_xy / 0.8);
      int nz = std::ceil(2 * max_z / 0.8);
      // Contract lattice by gamma factor in case of smearing where
      // smearing length is bound to the lattice cell length
      if (parameters_.smearing_mode == SmearingMode::Discrete ||
          parameters_.smearing_mode == SmearingMode::Triangular) {
        nz = static_cast<int>(std::ceil(2 * max_z / 0.8 * gam));
      }
      n_default = {n_xy, n_xy, nz};
    } else if (modus_.is_box()) {
      periodic_default = true;
      origin_default = {0., 0., 0.};
      const double bl = modus_.length();
      l_default = {bl, bl, bl};
      const int n_xyz = std::ceil(bl / 0.5);
      n_default = {n_xyz, n_xyz, n_xyz};
    } else if (modus_.is_sphere()) {
      // Maximal distance from (0, 0, 0) at which a particle
      // may be found at the end of the simulation
      const double max_d = modus_.radius() + end_time_;
      origin_default = {-max_d, -max_d, -max_d};
      l_default = {2 * max_d, 2 * max_d, 2 * max_d};
      // Go for approximately 0.8 fm cell size
      const int n_xyz = std::ceil(2 * max_d / 0.8);
      n_default = {n_xyz, n_xyz, n_xyz};
    }
    // Take lattice properties from config to assign them to all lattices
    const std::array<double, 3> l =
        config.take({"Lattice", "Sizes"}, l_default);
    const std::array<int, 3> n =
        config.take({"Lattice", "Cell_Number"}, n_default);
    const std::array<double, 3> origin =
        config.take({"Lattice", "Origin"}, origin_default);
    const bool periodic =
        config.take({"Lattice", "Periodic"}, periodic_default);

    logg[LExperiment].info()
        << "Lattice is ON. Origin = (" << origin[0] << "," << origin[1] << ","
        << origin[2] << "), sizes = (" << l[0] << "," << l[1] << "," << l[2]
        << "), number of cells = (" << n[0] << "," << n[1] << "," << n[2]
        << ")";

    if (printout_lattice_td_ || printout_full_lattice_any_td_) {
      dens_type_lattice_printout_ = output_parameters.td_dens_type;
      printout_tmn_ = output_parameters.td_tmn;
      printout_tmn_landau_ = output_parameters.td_tmn_landau;
      printout_v_landau_ = output_parameters.td_v_landau;
      printout_j_QBS_ = output_parameters.td_jQBS;
    }
    if (printout_tmn_ || printout_tmn_landau_ || printout_v_landau_) {
      Tmn_ = make_unique<RectangularLattice<EnergyMomentumTensor>>(
          l, n, origin, periodic, LatticeUpdate::AtOutput);
    }
    if (printout_j_QBS_) {
      j_QBS_lat_ = make_unique<DensityLattice>(l, n, origin, periodic,
                                               LatticeUpdate::AtOutput);
    }
    /* Create baryon and isospin density lattices regardless of config
       if potentials are on. This is because they allow to compute
       potentials faster */
    if (potentials_) {
      // Create auxiliary lattices for baryon four-current calculation
      old_jmu_auxiliary_ = make_unique<RectangularLattice<FourVector>>(
          l, n, origin, periodic, LatticeUpdate::EveryTimestep);
      new_jmu_auxiliary_ = make_unique<RectangularLattice<FourVector>>(
          l, n, origin, periodic, LatticeUpdate::EveryTimestep);
      four_gradient_auxiliary_ =
          make_unique<RectangularLattice<std::array<FourVector, 4>>>(
              l, n, origin, periodic, LatticeUpdate::EveryTimestep);

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
      if (potentials_->use_coulomb()) {
        jmu_el_lat_ = make_unique<DensityLattice>(l, n, origin, periodic,
                                                  LatticeUpdate::EveryTimestep);
        EM_lat_ = make_unique<
            RectangularLattice<std::pair<ThreeVector, ThreeVector>>>(
            l, n, origin, periodic, LatticeUpdate::EveryTimestep);
      }
      if (potentials_->use_vdf()) {
        jmu_B_lat_ = make_unique<DensityLattice>(l, n, origin, periodic,
                                                 LatticeUpdate::EveryTimestep);
        UB_lat_ = make_unique<RectangularLattice<FourVector>>(
            l, n, origin, periodic, LatticeUpdate::EveryTimestep);
        FB_lat_ = make_unique<
            RectangularLattice<std::pair<ThreeVector, ThreeVector>>>(
            l, n, origin, periodic, LatticeUpdate::EveryTimestep);
      }
      if (parameters_.field_derivatives_mode == FieldDerivativesMode::Direct) {
        // Create auxiliary lattices for field calculation
        old_fields_auxiliary_ = make_unique<RectangularLattice<FourVector>>(
            l, n, origin, periodic, LatticeUpdate::EveryTimestep);
        new_fields_auxiliary_ = make_unique<RectangularLattice<FourVector>>(
            l, n, origin, periodic, LatticeUpdate::EveryTimestep);
        fields_four_gradient_auxiliary_ =
            make_unique<RectangularLattice<std::array<FourVector, 4>>>(
                l, n, origin, periodic, LatticeUpdate::EveryTimestep);

        // Create the fields lattice
        fields_lat_ = make_unique<FieldsLattice>(l, n, origin, periodic,
                                                 LatticeUpdate::EveryTimestep);
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
  } else if (printout_lattice_td_ || printout_full_lattice_any_td_) {
    logg[LExperiment].error(
        "If you want Therm. VTK or Lattice output, configure a lattice for "
        "it.");
  } else if (potentials_ && potentials_->use_coulomb()) {
    logg[LExperiment].error(
        "Coulomb potential requires a lattice. Please add one to the "
        "configuration");
  }

  // Warning for the mean field calculation if lattice is not on.
  if ((potentials_ != nullptr) && (jmu_B_lat_ == nullptr)) {
    logg[LExperiment].warn() << "Lattice is NOT used. Mean-field energy is "
                             << "not going to be calculated.";
  }

  // Store pointers to potential and lattice accessible for Action
  if (parameters_.potential_affect_threshold) {
    UB_lat_pointer = UB_lat_.get();
    UI3_lat_pointer = UI3_lat_.get();
    pot_pointer = potentials_.get();
  }

  // Throw fatal if DerivativesMode == FiniteDifference and lattice is not on.
  if ((parameters_.derivatives_mode == DerivativesMode::FiniteDifference) &&
      (jmu_B_lat_ == nullptr)) {
    throw std::runtime_error(
        "Lattice is necessary to calculate finite difference gradients.");
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
 * Generate a string which will be printed to the screen when SMASH is running
 *
 * \param[in] ensembles The simulated particles: one Particles object per
 *            ensemble. The information about particles is used to check the
 *            conservation of the total energy and momentum as well as print
 *            other useful information.
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
std::string format_measurements(const std::vector<Particles> &ensembles,
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
 * \param[in] em_lattice Lattice containing the electric and magnetic field in
 * fm^-2 \param[in] parameters Parameters of the experiment, needed for the
 * access to the number of testparticles. \return Total mean field energy in the
 * Box.
 */
double calculate_mean_field_energy(
    const Potentials &potentials,
    RectangularLattice<smash::DensityOnLattice> &jmu_B_lat,
    RectangularLattice<std::pair<ThreeVector, ThreeVector>> *em_lattice,
    const ExperimentParameters &parameters);

/**
 * Generate the EventInfo object which is passed to outputs_.
 *
 * \param[in] ensembles The simulated particles: one Particles object per
 *            ensemble. Information about all particles (positions, momenta,
 *            etc.)is passed to the output.
 * \param[in] E_mean_field Value of the mean-field contribution to the total
 *            energy of the system at the current time.
 * \param[in] modus_impact_parameter The impact parameter
 * \param[in] parameters structure that holds various global parameters
 *            such as testparticle number, see \ref ExperimentParameters
 * \param[in] projectile_target_interact true if there was at least one
 *            collision
 * \param[in] kinematic_cut_for_SMASH_IC true if kinematic cuts in y or pT are
              enabled when exracting initial conditions for hydrodynamics
 */
EventInfo fill_event_info(const std::vector<Particles> &ensembles,
                          double E_mean_field, double modus_impact_parameter,
                          const ExperimentParameters &parameters,
                          bool projectile_target_interact,
                          bool kinematic_cut_for_SMASH_IC);

template <typename Modus>
void Experiment<Modus>::initialize_new_event() {
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

  for (Particles &particles : ensembles_) {
    particles.reset();
  }

  // Sample particles according to the initial conditions
  double start_time = -1.0;

  // Sample impact parameter only once per all ensembles
  // It should be the same for all ensembles
  if (modus_.is_collider()) {
    modus_.sample_impact();
    logg[LExperiment].info("Impact parameter = ", modus_.impact_parameter(),
                           " fm");
  }
  for (Particles &particles : ensembles_) {
    start_time = modus_.initial_conditions(&particles, parameters_);
  }
  /* For box modus make sure that particles are in the box. In principle, after
   * a correct initialization they should be, so this is just playing it safe.
   */
  for (Particles &particles : ensembles_) {
    modus_.impose_boundary_conditions(&particles, outputs_);
  }
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
  conserved_initial_ = QuantumNumbers(ensembles_);
  wall_actions_total_ = 0;
  previous_wall_actions_total_ = 0;
  interactions_total_ = 0;
  previous_interactions_total_ = 0;
  discarded_interactions_total_ = 0;
  total_pauli_blocked_ = 0;
  projectile_target_interact_.assign(parameters_.n_ensembles, false);
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
    // update_potentials();
    // if (parameters.outputclock->current_time() == 0.0 )
    // using the lattice is necessary
    if ((jmu_B_lat_ != nullptr)) {
      update_lattice(jmu_B_lat_.get(), old_jmu_auxiliary_.get(),
                     new_jmu_auxiliary_.get(), four_gradient_auxiliary_.get(),
                     LatticeUpdate::EveryTimestep, DensityType::Baryon,
                     density_param_, ensembles_,
                     parameters_.labclock->timestep_duration(), true);
      // Because there was no lattice at t=-Delta_t, the time derivatives
      // drho_dt and dj^mu/dt at t=0 are huge, while they shouldn't be; we
      // overwrite the time derivative to zero by hand.
      for (auto &node : *jmu_B_lat_) {
        node.overwrite_drho_dt_to_zero();
        node.overwrite_djmu_dt_to_zero();
      }
      E_mean_field = calculate_mean_field_energy(*potentials_, *jmu_B_lat_,
                                                 EM_lat_.get(), parameters_);
    }
  }
  initial_mean_field_energy_ = E_mean_field;
  logg[LExperiment].info() << format_measurements(
      ensembles_, 0u, conserved_initial_, time_start_,
      parameters_.labclock->current_time(), E_mean_field,
      initial_mean_field_energy_);

  // Output at event start
  for (const auto &output : outputs_) {
    for (int i_ens = 0; i_ens < parameters_.n_ensembles; i_ens++) {
      auto event_info = fill_event_info(
          ensembles_, E_mean_field, modus_.impact_parameter(), parameters_,
          projectile_target_interact_[i_ens], kinematic_cuts_for_IC_output_);
      output->at_eventstart(ensembles_[i_ens],
                            // Pretend each ensemble is an independent event
                            event_ * parameters_.n_ensembles + i_ens,
                            event_info);
    }
    // For thermodynamic output
    output->at_eventstart(ensembles_, event_);
    // For thermodynamic lattice output
    if (printout_full_lattice_any_td_) {
      switch (dens_type_lattice_printout_) {
        case DensityType::Baryon:
          output->at_eventstart(event_, ThermodynamicQuantity::EckartDensity,
                                DensityType::Baryon, *jmu_B_lat_);
          break;
        case DensityType::BaryonicIsospin:
          output->at_eventstart(event_, ThermodynamicQuantity::EckartDensity,
                                DensityType::BaryonicIsospin, *jmu_I3_lat_);
          break;
        case DensityType::None:
          break;
        default:
          output->at_eventstart(event_, ThermodynamicQuantity::EckartDensity,
                                DensityType::BaryonicIsospin, *jmu_custom_lat_);
      }
      if (printout_tmn_) {
        output->at_eventstart(event_, ThermodynamicQuantity::Tmn,
                              dens_type_lattice_printout_, *Tmn_);
      }
      if (printout_tmn_landau_) {
        output->at_eventstart(event_, ThermodynamicQuantity::TmnLandau,
                              dens_type_lattice_printout_, *Tmn_);
      }
      if (printout_v_landau_) {
        output->at_eventstart(event_, ThermodynamicQuantity::LandauVelocity,
                              dens_type_lattice_printout_, *Tmn_);
      }
      if (printout_j_QBS_) {
        output->at_eventstart(event_, ThermodynamicQuantity::j_QBS,
                              dens_type_lattice_printout_, *j_QBS_lat_);
      }
    }
  }

  /* In the ColliderModus, if Fermi motion is frozen, assign the beam momenta
   * to the nucleons in both the projectile and the target. Every ensemble
   * gets the same beam momenta, so no need to create beam_momenta_ vector
   * for every ensemble.
   */
  if (modus_.is_collider() && modus_.fermi_motion() == FermiMotion::Frozen) {
    for (ParticleData &particle : ensembles_[0]) {
      const double m = particle.effective_mass();
      double v_beam = 0.0;
      if (particle.belongs_to() == BelongsTo::Projectile) {
        v_beam = modus_.velocity_projectile();
      } else if (particle.belongs_to() == BelongsTo::Target) {
        v_beam = modus_.velocity_target();
      }
      const double gamma = 1.0 / std::sqrt(1.0 - v_beam * v_beam);
      beam_momentum_.emplace_back(
          FourVector(gamma * m, 0.0, 0.0, gamma * v_beam * m));
    }  // loop over particles
  }
}

template <typename Modus>
bool Experiment<Modus>::perform_action(Action &action, int i_ensemble,
                                       bool include_pauli_blocking) {
  Particles &particles = ensembles_[i_ensemble];
  // Make sure to skip invalid and Pauli-blocked actions.
  if (!action.is_valid(particles)) {
    discarded_interactions_total_++;
    logg[LExperiment].debug(~einhard::DRed(), "✘ ", action,
                           " (discarded: invalid)");
    return false;
  }
  action.generate_final_state();
  logg[LExperiment].debug("Process Type is: ", action.get_type());
  if (include_pauli_blocking && pauli_blocker_ &&
      action.is_pauli_blocked(ensembles_, *pauli_blocker_)) {
    total_pauli_blocked_++;
    return false;
  }

  // Prepare projectile_target_interact_, it's used for output
  // to signal that there was some interaction in this event
  if (modus_.is_collider()) {
    int count_target = 0, count_projectile = 0;
    for (const auto &incoming : action.incoming_particles()) {
      if (incoming.belongs_to() == BelongsTo::Projectile) {
        count_projectile++;
      } else if (incoming.belongs_to() == BelongsTo::Target) {
        count_target++;
      }
    }
    if (count_target > 0 && count_projectile > 0) {
      projectile_target_interact_[i_ensemble] = true;
    }
  }

  /* Make sure to pick a non-zero integer, because 0 is reserved for "no
   * interaction yet". */
  const auto id_process = static_cast<uint32_t>(interactions_total_ + 1);
  action.perform(&particles, id_process);
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
    // todo(oliiny): it's a rough density estimate from a single ensemble.
    // It might actually be appropriate for output. Discuss.
    rho = std::get<0>(current_eckart(r_interaction.threevec(), particles,
                                     density_param_, dens_type_, compute_grad,
                                     smearing));
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
void Experiment<Modus>::run_time_evolution(const double t_end,
                                           ActionList extra_actions) {
  while (parameters_.labclock->current_time() < t_end) {
    const double t = parameters_.labclock->current_time();
    const double dt =
        std::min(parameters_.labclock->timestep_duration(), t_end - t);
    logg[LExperiment].debug("Timestepless propagation for next ", dt, " fm/c.");


    ////////////////////////////////////////////////////////////////////////////////////////
    // Add demo action
    constexpr double r_x = 0.1;
    const FourVector pos_a{t, -r_x, 0., 0.};
    const FourVector pos_b{t, r_x, 0., 0.};

    ParticleData a{ParticleType::find(0x111)};  // pi0
    a.set_4position(pos_a);
    a.set_4momentum(a.pole_mass(), 1.0, 0., 0.);

    ParticleData b{ParticleType::find(0x111)};  // pi0
    b.set_4position(pos_b);
    b.set_4momentum(b.pole_mass(), -1.0, 0., 0.);

    ParticleList in_list{};
    ParticleList out_list{a, b};

    const double demo_act_time = t + 0.5 * dt;
    std::cout << "demo_act_time: " << demo_act_time << "\n";
    auto demo_act = make_unique<FreeforallAction>(in_list, out_list, demo_act_time);

    std::cout << demo_act->get_type() << '\n';

    // FreeforallAction demo_act(in_list, out_list, demo_act_time);



    // New idea directly perform action myself
    perform_action(*demo_act, 0);

    // std::vector<ActionPtr> extra_acts;
    // extra_acts.emplace_back(std::move(demo_act));
    //
    // actions[0].insert(std::move(extra_acts));
    // TODO only adds to first ensemble for now

    ////////////////////////////////////////////////////////////////////////////////////////

    // Perform forced thermalization if required
    if (thermalizer_ &&
        thermalizer_->is_time_to_thermalize(parameters_.labclock)) {
      const bool ignore_cells_under_treshold = true;
      // Thermodynamics in thermalizer is computed from all ensembles,
      // but thermalization actions act on each ensemble independently
      thermalizer_->update_thermalizer_lattice(ensembles_, density_param_,
                                               ignore_cells_under_treshold);
      const double current_t = parameters_.labclock->current_time();
      for (int i_ens = 0; i_ens < parameters_.n_ensembles; i_ens++) {
        thermalizer_->thermalize(ensembles_[i_ens], current_t,
                                 parameters_.testparticles);
        ThermalizationAction th_act(*thermalizer_, current_t);
        if (th_act.any_particles_thermalized()) {
          perform_action(th_act, i_ens);
        }
      }
    }

    std::vector<Actions> actions(parameters_.n_ensembles);
    for (int i_ens = 0; i_ens < parameters_.n_ensembles; i_ens++) {
      actions[i_ens].clear();
      if (ensembles_[i_ens].size() > 0 && action_finders_.size() > 0) {
        /* (1.a) Create grid. */
        const double min_cell_length = compute_min_cell_length(dt);
        logg[LExperiment].debug("Creating grid with minimal cell length ",
                                min_cell_length);
        /* For the hyper-surface-crossing actions also unformed particles are
         * searched and therefore needed on the grid. */
        const bool include_unformed_particles = IC_output_switch_;
        const auto &grid =
            use_grid_ ? modus_.create_grid(ensembles_[i_ens], min_cell_length,
                                           dt, parameters_.coll_crit,
                                           include_unformed_particles)
                      : modus_.create_grid(ensembles_[i_ens], min_cell_length,
                                           dt, parameters_.coll_crit,
                                           include_unformed_particles,
                                           CellSizeStrategy::Largest);

        const double gcell_vol = grid.cell_volume();
        /* (1.b) Iterate over cells and find actions. */
        grid.iterate_cells(
            [&](const ParticleList &search_list) {
              for (const auto &finder : action_finders_) {
                actions[i_ens].insert(finder->find_actions_in_cell(
                    search_list, dt, gcell_vol, beam_momentum_));
              }
            },
            [&](const ParticleList &search_list,
                const ParticleList &neighbors_list) {
              for (const auto &finder : action_finders_) {
                actions[i_ens].insert(finder->find_actions_with_neighbors(
                    search_list, neighbors_list, dt, beam_momentum_));
              }
            });
      }
    }

    /* \todo (optimizations) Adapt timestep size here */

    /* (2) Propagate from action to action until next output or timestep end */
    const double end_timestep_time =
        std::min(parameters_.labclock->next_time(), t_end);
    while (next_output_time() <= end_timestep_time) {
      for (int i_ens = 0; i_ens < parameters_.n_ensembles; i_ens++) {
        run_time_evolution_timestepless(actions[i_ens], i_ens,
                                        next_output_time(), t_end);
      }
      ++(*parameters_.outputclock);

      // Avoid duplication of final output
      if (parameters_.outputclock->current_time() < t_end) {
        intermediate_output();
      }
    }
    for (int i_ens = 0; i_ens < parameters_.n_ensembles; i_ens++) {
      std::cout << "run_time_evolution_timestepless"
                << "\n";
      std::cout << "#ensemble: " << i_ens << "\n";
      std::cout << "actions[i_ens].size(): " << actions[i_ens].size() << "\n";
      run_time_evolution_timestepless(actions[i_ens], i_ens, end_timestep_time,
                                      t_end);
    }

    /* (3) Update potentials (if computed on the lattice) and
     *     compute new momenta according to equations of motion */
    if (potentials_) {
      update_potentials();
      update_momenta(ensembles_, parameters_.labclock->timestep_duration(),
                     *potentials_, FB_lat_.get(), FI3_lat_.get(),
                     EM_lat_.get());
    }

    /* (4) Expand universe if non-minkowskian metric; updates
     *     positions and momenta according to the selected expansion */
    if (metric_.mode_ != ExpansionMode::NoExpansion) {
      for (Particles &particles : ensembles_) {
        expand_space_time(&particles, parameters_, metric_);
      }
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
      std::string err_msg = conserved_initial_.report_deviations(ensembles_);
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
void Experiment<Modus>::propagate_and_shine(double to_time,
                                            Particles &particles) {
  const double dt =
      propagate_straight_line(&particles, to_time, beam_momentum_);
  if (dilepton_finder_ != nullptr) {
    for (const auto &output : outputs_) {
      dilepton_finder_->shine(particles, output.get(), dt);
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
void Experiment<Modus>::run_time_evolution_timestepless(
    Actions &actions, int i_ensemble, const double end_time_propagation,
    const double end_time_run) {
  Particles &particles = ensembles_[i_ensemble];
  logg[LExperiment].info("Timestepless propagation: ", "Actions size = ",
                         actions.size(), ", end time = ", end_time_propagation);

  // iterate over all actions
  while (!actions.is_empty()) {
    if (actions.earliest_time() > end_time_propagation) {
      logg[LExperiment].info(
          "actions.earliest_time() > end_time_propagation is true",
          actions.earliest_time(), end_time_propagation);
      break;
    }
    // get next action
    ActionPtr act = actions.pop();
    if (!act->is_valid(particles)) {
      discarded_interactions_total_++;
      logg[LExperiment].info(~einhard::DRed(), "✘ ", act,
                             " (discarded: invalid)");
      continue;
    }
    logg[LExperiment].info(~einhard::Green(), "✔ ", act,
                           ", action time = ", act->time_of_execution());

    /* (1) Propagate to the next action. */
    propagate_and_shine(act->time_of_execution(), particles);

    /* (2) Perform action.
     *
     * Update the positions of the incoming particles, because the information
     * in the action object will be outdated as the particles have been
     * propagated since the construction of the action. */
    act->update_incoming(particles);
    const bool performed = perform_action(*act, i_ensemble);

    /* No need to update actions for outgoing particles
     * if the action is not performed. */
    if (!performed) {
      continue;
    }

    /* (3) Update actions for newly-produced particles. */

    const double end_time_timestep =
        std::min(parameters_.labclock->next_time(), end_time_run);
    assert(!(end_time_propagation > end_time_timestep));
    // New actions are always search until the end of the current timestep
    const double time_left = end_time_timestep - act->time_of_execution();
    const ParticleList &outgoing_particles = act->outgoing_particles();
    // Grid cell volume set to zero, since there is no grid
    const double gcell_vol = 0.0;
    for (const auto &finder : action_finders_) {
      // Outgoing particles can still decay, cross walls...
      actions.insert(finder->find_actions_in_cell(outgoing_particles, time_left,
                                                  gcell_vol, beam_momentum_));
      // ... and collide with other particles.
      actions.insert(finder->find_actions_with_surrounding_particles(
          outgoing_particles, particles, time_left, beam_momentum_));
    }

    check_interactions_total(interactions_total_);
  }

  propagate_and_shine(end_time_propagation, particles);
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
  /// Auxiliary variable to communicate the time in the computational frame
  /// at the functions printing the thermodynamics lattice output
  double computational_frame_time = 0.0;
  if (potentials_) {
    // using the lattice is necessary
    if ((jmu_B_lat_ != nullptr)) {
      E_mean_field = calculate_mean_field_energy(*potentials_, *jmu_B_lat_,
                                                 EM_lat_.get(), parameters_);
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
      ensembles_, interactions_this_interval, conserved_initial_, time_start_,
      parameters_.outputclock->current_time(), E_mean_field,
      initial_mean_field_energy_);
  const LatticeUpdate lat_upd = LatticeUpdate::AtOutput;

  // save evolution data
  if (!(modus_.is_box() && parameters_.outputclock->current_time() <
                               modus_.equilibration_time())) {
    for (const auto &output : outputs_) {
      if (output->is_dilepton_output() || output->is_photon_output() ||
          output->is_IC_output()) {
        continue;
      }
      for (int i_ens = 0; i_ens < parameters_.n_ensembles; i_ens++) {
        auto event_info = fill_event_info(
            ensembles_, E_mean_field, modus_.impact_parameter(), parameters_,
            projectile_target_interact_[i_ens], kinematic_cuts_for_IC_output_);

        output->at_intermediate_time(ensembles_[i_ens], parameters_.outputclock,
                                     density_param_, event_info);
        computational_frame_time = event_info.current_time;
      }
      // For thermodynamic output
      output->at_intermediate_time(ensembles_, parameters_.outputclock,
                                   density_param_);

      // Thermodynamic output on the lattice versus time
      switch (dens_type_lattice_printout_) {
        case DensityType::Baryon:
          update_lattice(jmu_B_lat_.get(), lat_upd, DensityType::Baryon,
                         density_param_, ensembles_, false);
          output->thermodynamics_output(ThermodynamicQuantity::EckartDensity,
                                        DensityType::Baryon, *jmu_B_lat_);
          output->thermodynamics_lattice_output(*jmu_B_lat_,
                                                computational_frame_time);
          break;
        case DensityType::BaryonicIsospin:
          update_lattice(jmu_I3_lat_.get(), lat_upd,
                         DensityType::BaryonicIsospin, density_param_,
                         ensembles_, false);
          output->thermodynamics_output(ThermodynamicQuantity::EckartDensity,
                                        DensityType::BaryonicIsospin,
                                        *jmu_I3_lat_);
          output->thermodynamics_lattice_output(*jmu_I3_lat_,
                                                computational_frame_time);
          break;
        case DensityType::None:
          break;
        default:
          update_lattice(jmu_custom_lat_.get(), lat_upd,
                         dens_type_lattice_printout_, density_param_,
                         ensembles_, false);
          output->thermodynamics_output(ThermodynamicQuantity::EckartDensity,
                                        dens_type_lattice_printout_,
                                        *jmu_custom_lat_);
          output->thermodynamics_lattice_output(*jmu_custom_lat_,
                                                computational_frame_time);
      }
      if (printout_tmn_ || printout_tmn_landau_ || printout_v_landau_) {
        update_lattice(Tmn_.get(), lat_upd, dens_type_lattice_printout_,
                       density_param_, ensembles_, false);
        if (printout_tmn_) {
          output->thermodynamics_output(ThermodynamicQuantity::Tmn,
                                        dens_type_lattice_printout_, *Tmn_);
          output->thermodynamics_lattice_output(
              ThermodynamicQuantity::Tmn, *Tmn_, computational_frame_time);
        }
        if (printout_tmn_landau_) {
          output->thermodynamics_output(ThermodynamicQuantity::TmnLandau,
                                        dens_type_lattice_printout_, *Tmn_);
          output->thermodynamics_lattice_output(
              ThermodynamicQuantity::TmnLandau, *Tmn_,
              computational_frame_time);
        }
        if (printout_v_landau_) {
          output->thermodynamics_output(ThermodynamicQuantity::LandauVelocity,
                                        dens_type_lattice_printout_, *Tmn_);
          output->thermodynamics_lattice_output(
              ThermodynamicQuantity::LandauVelocity, *Tmn_,
              computational_frame_time);
        }
      }
      if (EM_lat_) {
        output->fields_output("Efield", "Bfield", *EM_lat_);
      }
      if (printout_j_QBS_) {
        output->thermodynamics_lattice_output(
            *j_QBS_lat_, computational_frame_time, ensembles_, density_param_);
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
      update_lattice(jmu_I3_lat_.get(), old_jmu_auxiliary_.get(),
                     new_jmu_auxiliary_.get(), four_gradient_auxiliary_.get(),
                     LatticeUpdate::EveryTimestep, DensityType::BaryonicIsospin,
                     density_param_, ensembles_,
                     parameters_.labclock->timestep_duration(), true);
    }
    if ((potentials_->use_skyrme() || potentials_->use_symmetry()) &&
        jmu_B_lat_ != nullptr) {
      update_lattice(jmu_B_lat_.get(), old_jmu_auxiliary_.get(),
                     new_jmu_auxiliary_.get(), four_gradient_auxiliary_.get(),
                     LatticeUpdate::EveryTimestep, DensityType::Baryon,
                     density_param_, ensembles_,
                     parameters_.labclock->timestep_duration(), true);
      const size_t UBlattice_size = UB_lat_->size();
      for (size_t i = 0; i < UBlattice_size; i++) {
        auto jB = (*jmu_B_lat_)[i];
        const FourVector flow_four_velocity_B =
            std::abs(jB.rho()) > very_small_double ? jB.jmu_net() / jB.rho()
                                                   : FourVector();
        double baryon_density = jB.rho();
        ThreeVector baryon_grad_j0 = jB.grad_j0();
        ThreeVector baryon_dvecj_dt = jB.dvecj_dt();
        ThreeVector baryon_curl_vecj = jB.curl_vecj();
        if (potentials_->use_skyrme()) {
          (*UB_lat_)[i] =
              flow_four_velocity_B * potentials_->skyrme_pot(baryon_density);
          (*FB_lat_)[i] =
              potentials_->skyrme_force(baryon_density, baryon_grad_j0,
                                        baryon_dvecj_dt, baryon_curl_vecj);
        }
        if (potentials_->use_symmetry() && jmu_I3_lat_ != nullptr) {
          auto jI3 = (*jmu_I3_lat_)[i];
          const FourVector flow_four_velocity_I3 =
              std::abs(jI3.rho()) > very_small_double
                  ? jI3.jmu_net() / jI3.rho()
                  : FourVector();
          (*UI3_lat_)[i] = flow_four_velocity_I3 *
                           potentials_->symmetry_pot(jI3.rho(), baryon_density);
          (*FI3_lat_)[i] = potentials_->symmetry_force(
              jI3.rho(), jI3.grad_j0(), jI3.dvecj_dt(), jI3.curl_vecj(),
              baryon_density, baryon_grad_j0, baryon_dvecj_dt,
              baryon_curl_vecj);
        }
      }
    }
    if (potentials_->use_coulomb()) {
      update_lattice(jmu_el_lat_.get(), LatticeUpdate::EveryTimestep,
                     DensityType::Charge, density_param_, ensembles_, true);
      for (size_t i = 0; i < EM_lat_->size(); i++) {
        ThreeVector electric_field = {0., 0., 0.};
        ThreeVector position = jmu_el_lat_->cell_center(i);
        jmu_el_lat_->integrate_volume(electric_field,
                                      Potentials::E_field_integrand,
                                      potentials_->coulomb_r_cut(), position);
        ThreeVector magnetic_field = {0., 0., 0.};
        jmu_el_lat_->integrate_volume(magnetic_field,
                                      Potentials::B_field_integrand,
                                      potentials_->coulomb_r_cut(), position);
        (*EM_lat_)[i] = std::make_pair(electric_field, magnetic_field);
      }
    }  // if ((potentials_->use_skyrme() || ...
    if (potentials_->use_vdf() && jmu_B_lat_ != nullptr) {
      update_lattice(jmu_B_lat_.get(), old_jmu_auxiliary_.get(),
                     new_jmu_auxiliary_.get(), four_gradient_auxiliary_.get(),
                     LatticeUpdate::EveryTimestep, DensityType::Baryon,
                     density_param_, ensembles_,
                     parameters_.labclock->timestep_duration(), true);
      if (parameters_.field_derivatives_mode == FieldDerivativesMode::Direct) {
        update_fields_lattice(
            fields_lat_.get(), old_fields_auxiliary_.get(),
            new_fields_auxiliary_.get(), fields_four_gradient_auxiliary_.get(),
            jmu_B_lat_.get(), LatticeUpdate::EveryTimestep, *potentials_,
            parameters_.labclock->timestep_duration());
      }
      const size_t UBlattice_size = UB_lat_->size();
      for (size_t i = 0; i < UBlattice_size; i++) {
        auto jB = (*jmu_B_lat_)[i];
        (*UB_lat_)[i] = potentials_->vdf_pot(jB.rho(), jB.jmu_net());
        switch (parameters_.field_derivatives_mode) {
          case FieldDerivativesMode::ChainRule:
            (*FB_lat_)[i] = potentials_->vdf_force(
                jB.rho(), jB.drho_dxnu().x0(), jB.drho_dxnu().threevec(),
                jB.grad_rho_cross_vecj(), jB.jmu_net().x0(), jB.grad_j0(),
                jB.jmu_net().threevec(), jB.dvecj_dt(), jB.curl_vecj());
            break;
          case FieldDerivativesMode::Direct:
            auto Amu = (*fields_lat_)[i];
            (*FB_lat_)[i] = potentials_->vdf_force(
                Amu.grad_A0(), Amu.dvecA_dt(), Amu.curl_vecA());
            break;
        }
      }  // for (size_t i = 0; i < UBlattice_size; i++)
    }    // if potentials_->use_vdf()
  }
}

template <typename Modus>
void Experiment<Modus>::do_final_decays() {
  /* At end of time evolution: Force all resonances to decay. In order to handle
   * decay chains, we need to loop until no further actions occur. */
  bool actions_performed, decays_found;
  uint64_t interactions_old;
  do {
    decays_found = false;
    interactions_old = interactions_total_;
    for (int i_ens = 0; i_ens < parameters_.n_ensembles; i_ens++) {
      Actions actions;

      // Dileptons: shining of remaining resonances
      if (dilepton_finder_ != nullptr) {
        for (const auto &output : outputs_) {
          dilepton_finder_->shine_final(ensembles_[i_ens], output.get(), true);
        }
      }
      // Find actions.
      for (const auto &finder : action_finders_) {
        auto found_actions = finder->find_final_actions(ensembles_[i_ens]);
        if (!found_actions.empty()) {
          actions.insert(std::move(found_actions));
          decays_found = true;
        }
      }
      // Perform actions.
      while (!actions.is_empty()) {
        perform_action(*actions.pop(), i_ens, false);
      }
    }
    actions_performed = interactions_total_ > interactions_old;
    // Throw an error if actions were found but not performed
    if (decays_found && !actions_performed) {
      throw std::runtime_error("Final decays were found but not performed.");
    }
    // loop until no more decays occur
  } while (actions_performed);

  // Dileptons: shining of stable particles at the end
  if (dilepton_finder_ != nullptr) {
    for (const auto &output : outputs_) {
      for (Particles &particles : ensembles_) {
        dilepton_finder_->shine_final(particles, output.get(), false);
      }
    }
  }
}

template <typename Modus>
void Experiment<Modus>::final_output() {
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
        E_mean_field = calculate_mean_field_energy(*potentials_, *jmu_B_lat_,
                                                   EM_lat_.get(), parameters_);
      }
    }
    logg[LExperiment].info() << format_measurements(
        ensembles_, interactions_this_interval, conserved_initial_, time_start_,
        end_time_, E_mean_field, initial_mean_field_energy_);
    int total_particles = 0;
    for (const Particles &particles : ensembles_) {
      total_particles += particles.size();
    }
    if (IC_output_switch_ && (total_particles == 0)) {
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
               "the assumption for the stochastic criterion of\n1 interaction "
               "per particle per timestep is probably violated. Consider "
               "reducing the timestep size.";
      }

      logg[LExperiment].info() << "Final interaction number: "
                               << interactions_total_ - wall_actions_total_;
    }

    // Check if there are unformed particles
    int unformed_particles_count = 0;
    for (const Particles &particles : ensembles_) {
      for (const ParticleData &particle : particles) {
        if (particle.formation_time() > end_time_) {
          unformed_particles_count++;
        }
      }
    }
    if (unformed_particles_count > 0) {
      logg[LExperiment].warn(
          "End time might be too small. ", unformed_particles_count,
          " unformed particles were found at the end of the evolution.");
    }
  }

  // Keep track of how many ensembles had interactions
  count_nonempty_ensembles();

  for (const auto &output : outputs_) {
    for (int i_ens = 0; i_ens < parameters_.n_ensembles; i_ens++) {
      auto event_info = fill_event_info(
          ensembles_, E_mean_field, modus_.impact_parameter(), parameters_,
          projectile_target_interact_[i_ens], kinematic_cuts_for_IC_output_);
      output->at_eventend(ensembles_[i_ens],
                          // Pretend each ensemble is an independent event
                          event_ * parameters_.n_ensembles + i_ens, event_info);
    }
    // For thermodynamic output
    output->at_eventend(ensembles_, event_);

    // For thermodynamic lattice output
    if (dens_type_lattice_printout_ != DensityType::None) {
      output->at_eventend(event_, ThermodynamicQuantity::EckartDensity,
                          dens_type_lattice_printout_);
      output->at_eventend(ThermodynamicQuantity::EckartDensity);
    }
    if (printout_tmn_) {
      output->at_eventend(event_, ThermodynamicQuantity::Tmn,
                          dens_type_lattice_printout_);
      output->at_eventend(ThermodynamicQuantity::Tmn);
    }
    if (printout_tmn_landau_) {
      output->at_eventend(event_, ThermodynamicQuantity::TmnLandau,
                          dens_type_lattice_printout_);
      output->at_eventend(ThermodynamicQuantity::TmnLandau);
    }
    if (printout_v_landau_) {
      output->at_eventend(event_, ThermodynamicQuantity::LandauVelocity,
                          dens_type_lattice_printout_);
      output->at_eventend(ThermodynamicQuantity::LandauVelocity);
    }
    if (printout_j_QBS_) {
      output->at_eventend(ThermodynamicQuantity::j_QBS);
    }
  }
}

template <typename Modus>
void Experiment<Modus>::count_nonempty_ensembles() {
  for (bool has_interaction : projectile_target_interact_) {
    if (has_interaction) {
      nonempty_ensembles_++;
    }
  }
}

template <typename Modus>
bool Experiment<Modus>::is_finished() {
  if (event_counting_ == EventCounting::FixedNumber) {
    return event_ >= nevents_;
  }
  if (event_counting_ == EventCounting::MinimumNonEmpty) {
    if (nonempty_ensembles_ >= minimum_nonempty_ensembles_) {
      return true;
    }
    if (event_ >= max_events_) {
      logg[LExperiment].warn()
          << "Maximum number of events (" << max_events_
          << ") exceeded. Stopping calculation. "
          << "The fraction of empty ensembles is "
          << (1.0 - static_cast<double>(nonempty_ensembles_) /
                        (event_ * parameters_.n_ensembles))
          << ". If this fraction is expected, try increasing the "
             "Maximum_Ensembles_Run.";
      return true;
    }
    return false;
  }
  throw std::runtime_error("Event counting option is invalid");
  return false;
}

template <typename Modus>
void Experiment<Modus>::run() {
  const auto &mainlog = logg[LMain];
  for (event_ = 0; !is_finished(); event_++) {
    mainlog.info() << "Event " << event_;

    // Sample initial particles, start clock, some printout and book-keeping
    initialize_new_event();

    run_time_evolution(end_time_, {});

    if (force_decays_) {
      do_final_decays();
    }

    // Output at event end
    final_output();
  }
}

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_EXPERIMENT_H_
