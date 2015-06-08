/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ACTION_H_
#define SRC_INCLUDE_ACTION_H_

#include <stdexcept>
#include <vector>

#include "pauliblocking.h"
#include "particles.h"
#include "processbranch.h"
#include "random.h"

namespace Smash {

/**
 * \ingroup action
 * Action is the base class for a generic process that takes a number of
 * incoming particles and transforms them into any number of outgoing particles.
 * Currently such an action can be either a decay or a two-body collision
 * (see derived classes).
 */
class Action {
 public:
  /**
   * Construct an action object.
   *
   * \param[in] in_part list of incoming particles
   * \param[in] time Time at which the action is supposed to take place
   */
  Action(const ParticleList &in_part, float time);

  /** Copying is disabled. Use pointers or create a new Action. */
  Action(const Action &) = delete;

  /** Virtual Destructor.
   * The declaration of the destructor is necessary to make it virtual.
   */
  virtual ~Action();

  /** For sorting by time of execution. */
  bool operator<(const Action &rhs) const {
    return time_of_execution_ < rhs.time_of_execution_;
  }

  /** Return the raw weight value, which is a cross section in case of a
   * ScatterAction and a decay width in case of a DecayAction.
   *
   * Prefer to use a more specific function.
   */
  virtual float raw_weight_value() const = 0;

  /** Return the process type. */
  ProcessType get_type() const {
    return process_type_;
  }

  /** Add a new subprocess.  */
  template<typename Branch>
  void add_process(ProcessBranchPtr<Branch> &p,
                   ProcessBranchList<Branch>& subprocesses,
                   float& total_weight) {
    if (p->weight() > really_small) {
      total_weight += p->weight();
      subprocesses.emplace_back(std::move(p));
    }
  }
  /** Add several new subprocesses at once.  */
  template<typename Branch>
  void add_processes(ProcessBranchList<Branch> pv,
      ProcessBranchList<Branch>& subprocesses, float& total_weight) {
    if (subprocesses.empty()) {
      subprocesses = std::move(pv);
      for (auto &proc : subprocesses) {
        total_weight += proc->weight();
      }
    } else {
      subprocesses.reserve(subprocesses.size() + pv.size());
      for (auto &proc : pv) {
        total_weight += proc->weight();
        subprocesses.emplace_back(std::move(proc));
      }
    }
  }

  /**
   * Generate the final state for this action.
   *
   * This function selects a subprocess by Monte-Carlo decision and sets up
   * the final-state particles in phase space.
   */
  virtual void generate_final_state() = 0;

  /**
   * Actually perform the action, e.g. carry out a decay or scattering by
   * updating the particle list.
   *
   * This function removes the initial-state particles from the particle list
   * and then inserts the final-state particles. It does not do any sanity
   * checks, but assumes that is_valid has been called to determine if the
   * action is still valid.
   */
  virtual void perform(Particles *particles, size_t &id_process);

  /**
   * Check whether the action still applies.
   *
   * It can happen that a different action removed the incoming_particles from
   * the set of existing particles in the experiment, or that the particle has
   * scattered elastically in the meantime. In this case the Action doesn't
   * apply anymore and should be discarded.
   */
  bool is_valid(const Particles &) const;

  /**
   *  Check if the action is Pauli-blocked. If there are baryons in the final
   *  state then blocking probability is \f$ 1 - \Pi (1-f_i) \f$, where
   *  product is taken by all fermions in the final state and \f$ f_i \f$
   *  denotes phase-space density at the position of i-th final-state
   *  fermion.
   */
  bool is_pauli_blocked(const Particles &, const PauliBlocker &) const;

  /**
   * Return the list of particles that go into the interaction.
   */
  const ParticleList& incoming_particles() const;

  /**
   * Return the list of particles that resulted from the interaction.
   */
  const ParticleList &outgoing_particles() const { return outgoing_particles_; }

  /**
   * Return the time at which the action is supposed to be performed.
   */
  float time_of_execution() const { return time_of_execution_; }

  /**
   * Convert the time of execution from relative time to global time.
   *
   * This function should only be called once per action.
   *
   * \param current_global_time The current global time
   */
  void make_time_global(float current_global_time) {
    time_of_execution_ += current_global_time;
  }

  /** Check various conservation laws. */
  void check_conservation(const size_t &id_process) const;

  /** Get the interaction point */
  FourVector get_interaction_point();

  /**
   * \ingroup exception
   * Thrown for example when ScatterAction is called to perform with a wrong
   * number of final-state particles or when the energy is too low to produce
   * the resonance.
   */
  class InvalidResonanceFormation : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };

 protected:
  /** List with data of incoming particles.  */
  ParticleList incoming_particles_;
  /**
   * Initially this stores only the PDG codes of final-state particles.
   *
   * After perform was called it contains the complete particle data of the
   * outgoing particles.
   */
  ParticleList outgoing_particles_;
  /** time at which the action is supposed to be performed  */
  float time_of_execution_;
  /** type of process */
  ProcessType process_type_;

  /// determine the total energy in the center-of-mass frame
  /// \fpPrecision Why \c double?
  virtual double sqrt_s() const = 0;

  /**
   * Decide for a particular final-state channel via Monte-Carlo
   * and return it as a ProcessBranch
   */
  template<typename Branch>
  const Branch* choose_channel(
      const ProcessBranchList<Branch>& subprocesses, float total_weight) {
    const auto &log = logger<LogArea::Action>();
    float random_weight = Random::uniform(0.f, total_weight);
    float weight_sum = 0.;
    /* Loop through all subprocesses and select one by Monte Carlo, based on
     * their weights.  */
    for (const auto &proc : subprocesses) {
      /* All processes apart from strings should have
       * a well-defined final state. */
      if (proc->particle_number() < 1
          && proc->get_type() != ProcessType::String) {
        continue;
      }
      weight_sum += proc->weight();
      if (random_weight <= weight_sum) {
        /* Return the full process information. */
         return proc.get();
      }
    }
    /* Should never get here. */
    log.fatal(source_location, "Problem in choose_channel: ",
              subprocesses.size(), " ", weight_sum, " ", total_weight, " ",
    //          random_weight, "\n", *this);
              random_weight, "\n");
    throw std::runtime_error("problem in choose_channel");
  }

  /**
   * Sample final state momenta (and masses) in general X->2 processes,
   * using an isotropical angular distribution.
   *
   * \throws InvalidResonanceFormation
   */
  virtual void sample_cms_momenta();

  /**
   * \ingroup logging
   * Writes information about this action to the \p out stream.
   */
  virtual void format_debug_output(std::ostream &out) const = 0;

  /**
   * \ingroup logging
   * Dispatches formatting to the virtual Action::format_debug_output function.
   */
  friend std::ostream &operator<<(std::ostream &out, const Action &action) {
    action.format_debug_output(out);
    return out;
  }
};


inline std::vector<ActionPtr> &operator+=(std::vector<ActionPtr> &lhs,
                                          std::vector<ActionPtr> &&rhs) {
  if (lhs.size() == 0) {
    lhs = std::move(rhs);
  } else {
    lhs.insert(lhs.end(), std::make_move_iterator(rhs.begin()),
               std::make_move_iterator(rhs.end()));
  }
  return lhs;
}

/**
 * \ingroup logging
 * Convenience: dereferences the ActionPtr to Action.
 */
inline std::ostream &operator<<(std::ostream &out, const ActionPtr &action) {
  return out << *action;
}

/**
 * \ingroup logging
 * Writes multiple actions to the \p out stream.
 */
std::ostream &operator<<(std::ostream &out, const ActionList &actions);

}  // namespace Smash

#endif  // SRC_INCLUDE_ACTION_H_
