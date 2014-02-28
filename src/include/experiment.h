/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_EXPERIMENT_H_
#define SRC_INCLUDE_EXPERIMENT_H_

#include <memory>
#include <list>
#include <stdexcept>
#include <string>

#include "include/crosssections.h"
#include "include/experimentparameters.h"
#include "include/modusdefault.h"
#include "include/outputroutines.h"
#include "include/parameters.h"
#include "include/particles.h"

#ifndef DOXYGEN
namespace boost {
namespace filesystem {
class path;
}  // namespace filesystem
}  // namespace boost
#endif

namespace Smash {
class Configuration;

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
  virtual ~ExperimentBase() {}

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
   */
  static std::unique_ptr<ExperimentBase> create(Configuration &config);

  virtual void run(const boost::filesystem::path &path) = 0;

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
  virtual void run(const boost::filesystem::path &path) override;

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

  void initialize(const boost::filesystem::path &path);
  void run_time_evolution();
  void end();

  void print_startup();

  float energy_total(Particles *particles);

  inline timespec set_timer_start();

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
   * Struct of several member variables.
   * These variables are combined into a struct for efficient input to functions
   * outside of this class.
   */
  ExperimentParameters parameters_;

  /**
   * ?
   *
   * \todo CrossSections needs a rename?
   */
  CrossSections cross_sections_;

  /**
   * Number of events.
   *
   * \todo Explain what event means
   * \todo What does the number of events imply for the experiment?
   */
  int nevents_ = 0;

  /// number of steps
  int steps_ = 10000;
  /// number of steps before giving measurables
  int output_interval_ = 100;
  /// initial seed_ for random generator
  /// \todo Why is it a member? It's read, then used for seeding and
  ///       then never needed again, no?
  int64_t seed_ = 1;
  /// initial total energy of the system
  float energy_initial_ = 0.f;
  /// starting time of the simulation
  timespec time_start_ = set_timer_start();
};

}  // namespace Smash

#endif  // SRC_INCLUDE_EXPERIMENT_H_
