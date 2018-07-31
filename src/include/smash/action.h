/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ACTION_H_
#define SRC_INCLUDE_ACTION_H_

#include <stdexcept>
#include <utility>
#include <vector>

#include "action_globals.h"
#include "lattice.h"
#include "particles.h"
#include "pauliblocking.h"
#include "potentials.h"
#include "processbranch.h"
#include "random.h"

namespace smash {

/**
 * \ingroup action
 * Action is the base class for a generic process that takes a number of
 * incoming particles and transforms them into any number of outgoing particles.
 * Currently such an action can be either a decay, a two-body collision, a
 * wallcrossing or a thermalization.
 * (see derived classes).
 */
class Action {
 public:
  /**
   * Construct an action object with incoming particles and relative time.
   *
   * \param[in] in_part list of incoming particles
   * \param[in] time time at which the action is supposed to take place
   *                 (relative to the current time of the incoming particles)
   */
  Action(const ParticleList &in_part, double time)
      : incoming_particles_(in_part),
        time_of_execution_(time + in_part[0].position().x0()) {}

  /**
   * Construct an action object with the incoming particles, relative time, and
   * the already known outgoing particles and type of the process.
   *
   * \param[in] in_part list of incoming particles
   * \param[in] out_part list of outgoing particles
   * \param[in] time time at which the action is supposed to take place
   *                 (relative to the current time of the incoming particles)
   * \param[in] type type of the interaction
   */
  Action(const ParticleData &in_part, const ParticleData &out_part, double time,
         ProcessType type)
      : incoming_particles_({in_part}),
        outgoing_particles_({out_part}),
        time_of_execution_(time + in_part.position().x0()),
        process_type_(type) {}

  /**
   * Construct an action object with the incoming particles, absolute time, and
   * the already known outgoing particles and type of the process.
   *
   * \param[in] in_part list of incoming particles
   * \param[in] out_part list of outgoing particles
   * \param[in] absolute_execution_time absolute time at which the action is
   *                                             supposed to take place
   * \param[in] type type of the interaction
   */
  Action(const ParticleList &in_part, const ParticleList &out_part,
         double absolute_execution_time, ProcessType type)
      : incoming_particles_(std::move(in_part)),
        outgoing_particles_(std::move(out_part)),
        time_of_execution_(absolute_execution_time),
        process_type_(type) {}

  /// Copying is disabled. Use pointers or create a new Action.
  Action(const Action &) = delete;

  /**
   * Virtual Destructor.
   * The declaration of the destructor is necessary to make it virtual.
   */
  virtual ~Action();

  /**
   * Determine whether one action takes place before another in time
   *
   * \return if the first argument action takes place before the other
   */
  bool operator<(const Action &rhs) const {
    return time_of_execution_ < rhs.time_of_execution_;
  }

  /**
   * Return the total weight value, which is mainly used for the weight
   * output entry. It has different meanings depending of the type of
   * action. It is the total cross section in case of a ScatterAction,
   * the total decay width in case of a DecayAction and the shining
   * weight in case of a DecayActionDilepton.
   *
   * Prefer to use a more specific function. If there is no weight for the
   * action type, 0 should be returned.
   *
   * \return total cross section, decay width or shining weight
   */
  virtual double get_total_weight() const = 0;

  /**
   * Return the specific weight for the chosen outgoing channel, which is mainly
   * used for the partial weight output entry. For scatterings it will be the
   * partial cross section, for decays (including dilepton decays) the partial
   * decay width.
   *
   * If there is no weight for the action type, 0 should be returned.
   *
   * \return specific weight for the chosen output channel.
   */
  virtual double get_partial_weight() const = 0;

  /**
   * Get the process type.
   *
   * \return type of the process
   */
  virtual ProcessType get_type() const { return process_type_; }

  /**
   * Add a new subprocess.
   *
   * \param[in] p process to be added
   * \param[out] subprocesses processes, where p is added to
   * \param[out] total_weight summed weights of all the subprocesses
   */
  template <typename Branch>
  void add_process(ProcessBranchPtr<Branch> &p,
                   ProcessBranchList<Branch> &subprocesses,
                   double &total_weight) {
    if (p->weight() > 0) {
      total_weight += p->weight();
      subprocesses.emplace_back(std::move(p));
    }
  }

  /**
   * Add several new subprocesses at once.
   *
   * \param[in] pv processes list to be added
   * \param[out] subprocesses processes, where pv are added to
   * \param[out] total_weight summed weights of all the subprocesses
   */
  template <typename Branch>
  void add_processes(ProcessBranchList<Branch> pv,
                     ProcessBranchList<Branch> &subprocesses,
                     double &total_weight) {
    subprocesses.reserve(subprocesses.size() + pv.size());
    for (auto &proc : pv) {
      if (proc->weight() > 0) {
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
   *
   * \param[in] id_process unique id of the performed process
   * \param[out] particles particle list that is updated
   *
   * Note that you are required to increase id_process before the next call,
   * such that you get unique numbers.
   */
  virtual void perform(Particles *particles, uint32_t id_process);

  /**
   * Check whether the action still applies.
   *
   * It can happen that a different action removed the incoming_particles from
   * the set of existing particles in the experiment, or that the particle has
   * scattered elastically in the meantime. In this case the Action doesn't
   * apply anymore and should be discarded.

   * \param[in] particles current particle list
   * \return true, if action still applies; false otherwise
   */
  bool is_valid(const Particles &particles) const;

  /**
   * Check if the action is Pauli-blocked.
   *
   * If there are baryons in the final
   * state then blocking probability is \f$ 1 - \Pi (1-f_i) \f$, where the
   * product is taken by all fermions in the final state and \f$ f_i \f$
   * denotes the phase-space density at the position of i-th final-state
   * fermion.
   *
   * \param[in] particles current particle list
   * \param[in] p_bl PauliBlocker that stores the configurations concerning
   *                                                              Pauli-blocking.
   * \return true, if the action is Pauli-blocked, false otherwise
   */
  bool is_pauli_blocked(const Particles &particles,
                        const PauliBlocker &p_bl) const;

  /**
   * Get the list of particles that go into the action.
   *
   * \return a list of incoming particles
   */
  const ParticleList &incoming_particles() const;

  /**
   * Update the incoming particles that are stored in this action to the state
   * they have in the global particle list.
   *
   * \param[in] particles current particle list
   */
  void update_incoming(const Particles &particles);

  /**
   * Get the list of particles that resulted from the action.
   *
   * \return list of outgoing particles
   */
  const ParticleList &outgoing_particles() const { return outgoing_particles_; }

  /**
   * Get the time at which the action is supposed to be performed
   *
   * \return absolute time in the calculation frame in fm/c
   */
  double time_of_execution() const { return time_of_execution_; }

  /**
   * Check various conservation laws.
   *
   * \param[in] id_process process id only used for debugging output
   */
  void check_conservation(const uint32_t id_process) const;

  /**
   * Determine the total energy in the center-of-mass frame [GeV]
   *
   * \return \f$ \sqrt{s}\f$ of incoming particles
   */
  double sqrt_s() const { return total_momentum().abs(); }

  /**
   * Calculate the total kinetic energy of the outgoing particles in
   * the center of mass frame in the presence (or absence) of the mean field
   * potentials.
   *
   * This function is used when the species of the outgoing
   * particles are already determined.
   *
   * \return total kinetic energy of the outgoing particles [GeV]
   */
  double kinetic_energy_cms() const;

  /**
   * Calculate the total kinetic energy of the outgoing particles in
   * the center of mass frame in the presence (or absence) of the mean field
   * potentials.
   *
   * This function is used to determine whether an action is
   * kinematically feasible.
   *
   * \tparam outs Type of outgoing particles
   * \param[in] potentials skyrme and asymmetry potential for particle [GeV]
   * \param[in] p_out_types outgoing particle types
   * \return total kinetic energy of the outgoing particles in
   *                                           the center of mass frame [GeV]
   */
  template <typename outs>
  double kinetic_energy_cms(std::pair<double, double> potentials,
                            outs p_out_types) const {
    /* scale_B returns the difference of the total force scales of the skyrme
     * potential between the initial and final states. */
    double scale_B = 0.0;
    /* scale_I3 returns the difference of the total force scales of the symmetry
     * potential between the initial and final states. */
    double scale_I3 = 0.0;
    for (const auto &p_in : incoming_particles_) {
      // Get the force scale of the incoming particle.
      const auto scale =
          ((pot_pointer != nullptr) ? pot_pointer->force_scale(p_in.type())
                                    : std::make_pair(0.0, 0));
      scale_B += scale.first;
      scale_I3 += scale.second * p_in.type().isospin3_rel();
    }
    for (const auto &p_out : p_out_types) {
      const auto scale = ((pot_pointer != nullptr)
                              ? pot_pointer->force_scale(type_of_pout(p_out))
                              : std::make_pair(0.0, 0));
      scale_B -= scale.first;
      scale_I3 -= scale.second * type_of_pout(p_out).isospin3_rel();
    }
    /* Rescale to get the potential difference between the
     * initial and final state.*/
    const double B_pot_diff = potentials.first * scale_B;
    const double I3_pot_diff = potentials.second * scale_I3;
    return sqrt_s() + B_pot_diff + I3_pot_diff;
  }

  /**
   * Get the interaction point
   *
   * \return four vector of interaction point
   */
  FourVector get_interaction_point() const;

  /**
   * Get the skyrme and asymmetry potential at the interaction point
   *
   * \return skyrme and asymmetry potential [GeV]
   */
  std::pair<double, double> get_potential_at_interaction_point() const;

  /**
   * Input the information on the potentials and store in global variables
   * defined in action_globals.
   *
   * \param[in] UB_lat skyrme potential on lattice
   * \param[in] UI3_lat symmmetry potential on lattice
   * \param[in] pot potential class
   */
  static void input_potential(RectangularLattice<double> *UB_lat,
                              RectangularLattice<double> *UI3_lat,
                              Potentials *pot);

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
  /// List with data of incoming particles.
  ParticleList incoming_particles_;

  /**
   * Initially this stores only the PDG codes of final-state particles.
   *
   * After perform was called it contains the complete particle data of the
   * outgoing particles.
   */
  ParticleList outgoing_particles_;

  /**
   * Time at which the action is supposed to be performed
   * (absolute time in the lab frame in fm/c).
   */
  const double time_of_execution_;

  /// type of process
  ProcessType process_type_;

  /// Sum of 4-momenta of incoming particles
  FourVector total_momentum() const {
    FourVector mom(0.0, 0.0, 0.0, 0.0);
    for (const auto &p : incoming_particles_) {
      mom += p.momentum();
    }
    return mom;
  }

  /**
   * Remove the sub-threshold processes from the list of sub processes.
   *
   * \tparam Branch Type of the processbranch
   * \param[out] subprocesses list of processes that are possible
   * \param[out] total_weight summed weight of all subprocess (after filtering)
   */
  template <typename Branch>
  void filter_channel(ProcessBranchList<Branch> &subprocesses,
                      double &total_weight) {
    const auto potentials = get_potential_at_interaction_point();
    /* Loop through all subprocesses and remove sub-threshold ones.*/
    for (auto proc = subprocesses.begin(); proc != subprocesses.end();) {
      /* Evaluate the total kinentic energy of the final state particles
       * of this new subprocess. */
      const auto out_particle_types = (*proc)->particle_types();
      const double kin_energy_cms = kinetic_energy_cms<ParticleTypePtrList>(
          potentials, out_particle_types);
      /* Reject the process if the total kinetic energy is smaller than the
       * threshold. */
      if (kin_energy_cms < (*proc)->threshold()) {
        total_weight -= (*proc)->weight();
        proc = subprocesses.erase(proc);
      } else {
        ++proc;
      }
    }
  }

  /**
   * Decide for a particular final-state channel via Monte-Carlo
   * and return it as a ProcessBranch

   * \tparam Branch Type of processbranch
   * \param[in] subprocesses list of possible processes
   * \param[in] total_weight summed weight of all processes
   * \return ProcessBranch that is sampled
   */
  template <typename Branch>
  const Branch *choose_channel(const ProcessBranchList<Branch> &subprocesses,
                               double total_weight) {
    const auto &log = logger<LogArea::Action>();
    double random_weight = random::uniform(0., total_weight);
    double weight_sum = 0.;
    /* Loop through all subprocesses and select one by Monte Carlo, based on
     * their weights.  */
    for (const auto &proc : subprocesses) {
      weight_sum += proc->weight();
      if (random_weight <= weight_sum) {
        /* Return the full process information. */
        return proc.get();
      }
    }
    /* Should never get here. */
    log.fatal(source_location,
              "Problem in choose_channel: ", subprocesses.size(), " ",
              weight_sum, " ", total_weight, " ",
              //          random_weight, "\n", *this);
              random_weight, "\n");
    abort();
    throw std::runtime_error("problem in choose_channel");
  }

  /**
   * Sample final-state masses in general X->2 processes
   * (thus also fixing the absolute c.o.m. momentum).
   *
   * \throws InvalidResonanceFormation
   * \return masses of final state particles
   */
  virtual std::pair<double, double> sample_masses() const;

  /**
   * Sample final-state angles in general X->2 processes
   * (here: using an isotropical angular distribution).
   *
   * \param[in] masses masses of each of the final state particles
   */
  virtual void sample_angles(std::pair<double, double> masses);

  /**
   * Sample the full 2-body phase-space (masses, momenta, angles)
   * in the center-of-mass frame for the final state particles.
   */
  void sample_2body_phasespace();

  /**
   * \ingroup logging
   * Writes information about this action to the \p out stream.
   *
   * \param[out] out out stream to be written to
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

 private:
  /**
   * Get the type of a given particle
   *
   * \param[in] p_out particle of which the type will be returned
   * \return type of given particle
   */
  const ParticleType &type_of_pout(const ParticleData &p_out) const {
    return p_out.type();
  }
  /**
   * Get the particle type for given pointer to a particle type.
   *
   * Helper function for kinetic_energy_cms
   *
   * \param[in] p_out pointer to a particle type
   * \return particle type
   */
  const ParticleType &type_of_pout(const ParticleTypePtr &p_out) const {
    return *p_out;
  }
};

/**
 * Append vector of action pointers
 *
 * \param[in] lhs vector of action pointers that is appended to
 * \param[in] rhs vector of action pointers that is appended
 * \return vector of action pointers containing lhs and rhs
 */
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

}  // namespace smash

#endif  // SRC_INCLUDE_ACTION_H_
