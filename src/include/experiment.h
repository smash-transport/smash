/*
 *    Copyright (c) 2013-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_EXPERIMENT_H_
#define SRC_INCLUDE_EXPERIMENT_H_

#include <boost/filesystem.hpp>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "actionfinderfactory.h"
#include "chrono.h"
#include "density.h"
#include "experimentparameters.h"
#include "forwarddeclarations.h"
#include "outputinterface.h"
#include "lattice.h"
#include "particles.h"
#include "pauliblocking.h"
#include "potentials.h"
#include "processbranch.h"
#include "quantumnumbers.h"
#include "random.h"

namespace Smash {

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
   * \param config The configuration object that sets all initial conditions of
   *               the experiment.
   * \param output_path The directory where the output files are written.
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
                                                bf::path output_path);

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
  void run() override;

 private:
  /**
   * Create a new Experiment.
   *
   * This constructor is only called from the ExperimentBase::create factory
   * method.
   *
   * \param config  The Configuration object contains all initial setup of the
   *                experiment. It is forwarded to the constructors of member
   *                variables as needed.
   *                Note that the object is passed by non-const reference. This
   *                is only necessary for bookkeeping: Values are not only read,
   *                but actually taken out of the object. Thus, all values that
   *                remain were not used.
   * \param output_path The directory where the output files are written.
   */
  explicit Experiment(Configuration config, bf::path output_path);

  /** Reads particle type information and cross sections information and
   * does the initialization of the system
   *
   * This is called in the beginning of each event.
   */
  void initialize_new_event();

  /** Perform the given action. */
  template <typename Container>
  void perform_action(const ActionPtr &action, size_t &interactions_total,
                      size_t &total_pauliblocked,
                      const Container &particles_before_actions);

  /** It generates the final state with the right kinematics and then writes
   * the given dilepton action in the dilepton output file, instead of
   * actually performing the action.
   */
  void write_dilepton_action(const ActionPtr &action,
                               const ParticleList &particles_before_actions);

  /** Runs the time evolution of an event
   *
   * Here, the time steps are looped over, collisions and decays are
   * carried out and particles are propagated.
   *
   * \param evt_num Running number of the event
   * \return The number of interactions from the event
   */
  size_t run_time_evolution(const int evt_num);

  /** Runs the time evolution of an event without time steps
   *
   * Here, all actions are looped over, collisions and decays are
   * carried out and particles are propagated.
   *
   * \param evt_num Running number of the event
   * \return The number of interactions from the event
   */
  size_t run_time_evolution_without_time_steps(const int evt_num);

  /** Performs the final decays of an event
   *
   * \param interactions_total The number of interactions so far
   */
  void do_final_decays(size_t &interactions_total);

  /** Output at the end of an event
   *
   * \param interactions_total The number of interactions from the event
   * \param evt_num Number of the event
   */
  void final_output(size_t interactions_total, const int evt_num);

  /** Intermediate output during an event
   *
   * \param evt_num Number of the event
   * \param interactions_total The total number of interactions so far
   * \param previous_interactions_total The number of interactions at the
   *                                    previous output
   */
  void intermediate_output(const int evt_num, size_t& interactions_total,
                           size_t& previous_interactions_total);

  /**
   * Struct of several member variables.
   * These variables are combined into a struct for efficient input to functions
   * outside of this class.
   */
  ExperimentParameters parameters_;

  /**
   * Instance of the Modus template parameter. May store modus-specific data
   * and contains modus-specific function implementations.
   */
  Modus modus_;

  /**
   * The particles interacting in the experiment.
   */
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
  std::unique_ptr<OutputInterface> dilepton_output_;

  /// The Action finder objects
  std::vector<std::unique_ptr<ActionFinderInterface>> action_finders_;

  /// The Dilepton Action Finder
  std::unique_ptr<ActionFinderInterface> dilepton_finder_;

  /// Lattices holding different physical quantities

  /** Baryon density, isospin projection density, custom density.
   *  In the config user asks for some kind of density for printout.
   *  Baryon and isospin projection density are anyway needed for potentials.
   *  If user asks for some other density type for printout, it will be handled
   *  using jmu_custom variable.
   */
  std::unique_ptr<DensityLattice> jmu_B_lat_, jmu_I3_lat_, jmu_custom_lat_;
  /// Type of density for lattice printout
  DensityType dens_type_lattice_printout_ = DensityType::none;

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
  const float end_time_;
  /** The clock's timestep size at start up
   *
   * Stored here so that the next event will remember this.
   */
  const float delta_time_startup_;

  /**
   * This indicates whether we force all resonances to decay in the last timestep.
   */
  const bool force_decays_;

  /**
   * This indicates whether to use the grid.
   */
  const bool use_grid_;

  /**
   * This indicates whether to use time steps.
   */
  const bool use_time_steps_;

  /** The conserved quantities of the system.
   *
   * This struct carries the sums of the single particle's various
   * quantities as measured at the beginning of the evolution and can be
   * used to regularly check if they are still good.
   *
   */
  QuantumNumbers conserved_initial_;
  /// system starting time of the simulation
  SystemTimePoint time_start_ = SystemClock::now();

  /// Type of density to be written to collision headers
  DensityType dens_type_;

  /**\ingroup logging
   * Writes the initial state for the Experiment to the output stream.
   * It automatically appends the output of the current Modus.
   */
  friend std::ostream &operator<< <>(std::ostream &out, const Experiment &e);
};
}  // namespace Smash

#endif  // SRC_INCLUDE_EXPERIMENT_H_
