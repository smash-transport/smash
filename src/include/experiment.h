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
#include "decayactionsfinderdilepton.h"
#include "energymomentumtensor.h"
#include "fourvector.h"
#include "grandcan_thermalizer.h"
#include "outputparameters.h"
#include "pauliblocking.h"
#include "potentials.h"
#include "propagation.h"
#include "quantumnumbers.h"
#include "thermalizationaction.h"

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
   * Reads particle type information and cross sections information and
   * does the initialization of the system
   *
   * This is called in the beginning of each event.
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

  /// Provides external access to particles_
  Particles* particles() { return &particles_; }

  /// Provides external access to modus_
  Modus* modus() { return &modus_; }

 private:
  /**
   * Perform the given action.
   *
   * \param[in] action The action to perform. If it performs, it'll modify
   *                   the private member particles_
   * \param[in] particles_before_actions A container with the ParticleData
   *                 from this time step before any actions were performed
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

  /// The particles interacting in the experiment.
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

  /// Lattices for Skyme potential gradients.
  std::unique_ptr<RectangularLattice<ThreeVector>> dUB_dr_lat_;

  /// Lattices for symmetry potential gradients.
  std::unique_ptr<RectangularLattice<ThreeVector>> dUI3_dr_lat_;

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
   *  Total number of interactions for current and for previous timestep.
   *  For timestepless mode the whole run time is considered as one timestep.
   */
  uint64_t interactions_total_ = 0, previous_interactions_total_ = 0;

  /**
   *  Total number of wall-crossings for current and for previous timestep.
   *  For timestepless mode the whole run time is considered as one timestep.
   */
  uint64_t wall_actions_total_ = 0, previous_wall_actions_total_ = 0;

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
}  // namespace smash

#endif  // SRC_INCLUDE_EXPERIMENT_H_
