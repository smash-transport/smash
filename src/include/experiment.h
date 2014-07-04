/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_EXPERIMENT_H_
#define SRC_INCLUDE_EXPERIMENT_H_

#include <memory>
#include <stdexcept>
#include <string>
#include <map>
#include <vector>

#include "chrono.h"
#include "crosssections.h"
#include "decayactionsfinder.h"
#include "experimentparameters.h"
#include "forwarddeclarations.h"
#include "outputinterface.h"
#include "particles.h"
#include "quantumnumbers.h"
#include "scatteractionsfinder.h"

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
   *
   * \return An owning pointer to the Experiment object, using the
   *         ExperimentBase interface.
   *
   * \throws InvalidModusRequest This exception is thrown if the \p
   *         Modus string in the \p config object does not contain a valid
   *         string.
   *
   * GENERAL:
   * --------
   */
  // !!USER:Input
  /**
   * \if user
   * \page input_general_ Input Section General
   * \endif
   *
   * `MODUS:` Choose a modus for the calculation, e.g.\ infinite matter
   * calculation, collision of two particles or collision of nuclei. The modus
   * will be configured in ref input_modi_. Recognized values are:
   *
   * \li `Nucleus` for collisions of nuclei or compound objects. See
   *     \if user
   *     \ref input_modi_nucleus_
   *     \else
   *     \ref NucleusModus
   *     \endif
   * \li `Sphere` for calculations of the expansion of a thermalized sphere.
   * See ref input_modi_sphere_
   * \li `Collider` ...
   * \li `Box` for infinite matter calculation in a rectangular box. See
   *     \if user
   *     \ref input_modi_box_
   *     \else
   *     \ref BoxModus
   *     \endif
   *
   * `DELTA_TIME:` Time step for the calculation, in fm/c.
   *
   * `END_TIME:` The time after which the evolution is stopped. Note
   * that the starting time depends on the chosen MODUS.
   *
   * 'OUTPUT:' List of the desired output format(s).
   *
   * `OUTPUT_INTERVAL:` Defines the period of intermediate output
   * of the status of the simulated system in Standard Output
   * and other output formats which support this functionality.
   *
   * `RANDOMSEED:` Initial seed for the random number generator. If this is
   * negative, the program starting time is used.
   *
   * `SIGMA:` Elastic cross-section.
   *
   * `TESTPARTICLES:` How many test particles per real particles should be
   * simulated.
   *
   * `NEVENTS:` Number of events to calculate.
   *
   * `particles:` ???
   *
   * `decaymodes:` ???
   */
  // !!/USER:Input
  static std::unique_ptr<ExperimentBase> create(Configuration &config);

  /**
   * Runs the experiment.
   *
   * The constructor does the setup of the experiment. The run function executes
   * the complete experiment.
   *
   * \param path The path where output files will be written to.
   */
  virtual void run(const bf::path &path) = 0;

  /**
   * Exception class that is thrown if an invalid modus is requested from the
   * Experiment factory.
   */
  struct InvalidModusRequest : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };
};

/**
 * The main class, where the simulation of an experiment is executed.
 *
 * The Experiment class is owns all data (maybe indirectly) relevant for the
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
  virtual void run(const bf::path &path) override;

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
   */
  explicit Experiment(Configuration &config);

  /** Reads particle type information and cross sections information and
   * does the initialization of the system
   *
   * This is called in the beginning of each event.
   */
  void initialize_new_event();

  /*
   *  Sets up SMASH outputs.
   */
  void set_outputs(const bf::path &path);

  /** Runs the time evolution of an event
   *
   * Here, the time steps are looped over, collisions and decays are
   * carried out and particles are propagated.
   *
   * \param evt_num Running number of the event
   */
  void run_time_evolution(const int evt_num);

  /** Output about the beginning of the event.
   *
   * This output goes to the console
   *
   * \param seed The random number engine seed
   */
  void print_startup(int64_t seed);

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
   * A list of output formaters. They will be called to write the state of the
   * particles to file.
   */
  OutputsList outputs_;

  /**
   * ?
   *
   * \todo CrossSections needs a rename?
   */
  CrossSections cross_sections_;

  /// The object that finds decays
  DecayActionsFinder decay_finder_;
  /// The object that finds scatterings
  ScatterActionsFinder scatter_finder_;

  /**
   * Number of events.
   *
   * \todo Explain what event means
   * \todo What does the number of events imply for the experiment?
   */
  const int nevents_;

  /// simulation time at which the evolution is stopped.
  const float end_time_ = 10.0f;
  /** The clock's timestep size at start up
   *
   * Stored here so that the next event will remember this.
   */
  const float delta_time_startup_;
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

  /** Output formats.
   *
   * Specifies the format of output; OSCAR, ROOT, VTK etc.
   */
  std::map<std::string, std::string> outputformats_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_EXPERIMENT_H_
