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
   * \param[in] time_of_execution time at which the action is supposed to take place
   */
  Action(const ParticleList &in_part, float time_of_execution);

  /// Copying is disabled. Use std::move or create a new Action.
  Action(const Action &) = delete;

  /**
   * Move constructor for Action.
   *
   * The move constructor makes moving efficient since it can move the three
   * std::vector member variables.
   */
  Action(Action &&);

  /** Virtual Destructor.
   * The declaration of the destructor is necessary to make it virtual.
   */
  virtual ~Action();

  /** For sorting by time of execution. */
  bool operator<(const Action &rhs) const {
    return time_of_execution_ < rhs.time_of_execution_;
  }

  /** Returns the total weight, which is a cross section in case of a ScatterAction
   * and a decay width in case of a DecayAction. */
  float weight() const {
    return total_weight_;
  }

  /** Return the process type. */
  ProcessBranch::ProcessType get_type() const {
    return process_type_;
  }

  /** Add a new subprocess.  */
  void add_process(ProcessBranch *p);
  /** Add several new subprocesses at once.  */
  void add_processes(ProcessBranchList pv);

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
  ParticleList incoming_particles() const;

  /**
   * Return the list of particles that resulted from the interaction.
   */
  const ParticleList &outgoing_particles() const { return outgoing_particles_; }

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
  /** list of possible subprocesses  */
  ProcessBranchList subprocesses_;
  /** time at which the action is supposed to be performed  */
  float time_of_execution_;
  /** sum of all subprocess weights  */
  float total_weight_;
  /** type of process */
  ProcessBranch::ProcessType process_type_;

  /// determine the total energy in the center-of-mass frame
  /// \fpPrecision Why \c double?
  virtual double sqrt_s() const = 0;

  /**
   * Decide for a particular final-state channel via Monte-Carlo
   * and return it as a ProcessBranch
   */
  const ProcessBranch* choose_channel();

  /**
   * Sample final state momenta (and masses) in general X->2 process.
   *
   * \throws InvalidResonanceFormation
   */
  void sample_cms_momenta();

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
