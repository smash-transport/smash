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
#include <string>

#include "include/crosssections.h"
#include "include/experimentparameters.h"
#include "include/modusdefault.h"
#include "include/outputroutines.h"
#include "include/parameters.h"
#include "include/particles.h"

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
  /* Virtual base class destructor
   * to avoid undefined behavior when destroying derived objects
   */
  virtual ~ExperimentBase() {}

  virtual void configure(std::list<Parameters> configuration) = 0;
  virtual void commandline_arg(int steps) = 0;

  static std::unique_ptr<ExperimentBase> create(std::string modus_chooser,
                                                int nevents);

  virtual void run(std::string path) = 0;
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
 public:
  explicit Experiment(int nevents)
      : particles_(nullptr), cross_sections_(nullptr), nevents_(nevents) {}

  virtual void configure(std::list<Parameters> configuration) override;
  virtual void commandline_arg(int steps) override;

  virtual void run(std::string path) override;

 private:
  void initialize(const char *path);
  void run_time_evolution();
  void end();

  void assign_params(std::list<Parameters> *configuration);

  void print_startup();

  float energy_total(Particles *particles);

  inline timespec set_timer_start();

  Modus modus_;
  Particles *particles_;
  CrossSections *cross_sections_;

  int nevents_;

  ExperimentParameters parameters_;

  /* number of steps */
  int steps_ = 10000;
  /* number of steps before giving measurables */
  int output_interval_ = 100;
  /* initial seed_ for random generator */
  int64_t seed_ = 1;
  /* initial total energy of the system */
  float energy_initial_ = 0.f;
  /* starting time of the simulation */
  timespec time_start_ = set_timer_start();
};

#endif  // SRC_INCLUDE_EXPERIMENT_H_
