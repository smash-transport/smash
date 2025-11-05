/*
 *
 *    Copyright (c) 2015-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_DECAYACTION_H_
#define SRC_INCLUDE_SMASH_DECAYACTION_H_

#include <optional>
#include <utility>

#include "action.h"

namespace smash {

/**
 * \ingroup action
 * DecayAction is a special action which takes one single particle in the
 * initial state and makes it decay into a number of daughter particles
 * (currently two or three).
 */
class DecayAction : public Action {
 public:
  /**
   * Construct a DecayAction from a particle \p p.
   *
   * It does not initialize the list of possible decay processes. You need to
   * call add_processes after construction.
   *
   * \param[in] p The particle that should decay if the action is performed.
   * \param[in] time Time at which the action is supposed to take place
   * \param[in] spin_interaction_type Which type of spin interaction to use
   */
  DecayAction(
      const ParticleData &p, double time,
      SpinInteractionType spin_interaction_type = SpinInteractionType::Off);

  /**
   * Add several new decays at once.
   * \param[in] pv List of decays to be added.
   */
  void add_decays(DecayBranchList pv);

  /**
   * Add one new decay.
   * \param[in] p Decay to be added.
   */
  void add_decay(DecayBranchPtr p);

  /**
   * Generate the final state of the decay process.
   * Performs a decay of one particle to two or three particles.
   *
   * \throws InvalidDecay
   */
  void generate_final_state() override;

  /**
   * Sample the masses of the final particles
   * \returns Pair of sampled masses of particle 1 and 2
   */
  std::pair<double, double> sample_masses(
      double kinetic_energy_cm) const override;

  /**
   * Sample the full 2-body phase space (masses, momenta, angles)
   * in the center-of-mass frame for the final-state particles.
   *
   * \see Action::sample_2body_phasespace for the base behaviour.
   *
   * This overrides the base implementation to integrate with DecayAction’s
   * channel selection and potential-aware kinematics. If sample_masses()
   * returns NaNs (i.e., the previously chosen channel is kinematically
   * forbidden once potentials are considered), this method sets
   * was_2body_phase_space_sampled_with_potentials_as_valid_ to 'false' so the
   * caller can fall back to another channel instead of throwing. Otherwise it
   * sets it to 'true' and proceeds with momentum/angle sampling (which may use
   * L_).
   */
  void sample_2body_phasespace() override;

  /// Return the total width of the decay process.
  double get_total_weight() const override { return total_width_; }

  /**
   * Get partial width of chosen channel
   * \return Partial width of chosen channel
   */
  double get_partial_weight() const override { return partial_width_; }

  /**
   * Get total decay width
   * \return Total width of decay
   */
  double total_width() const { return total_width_; }

  /**
   * \ingroup exception
   * Thrown when DecayAction is called to perform with 0 or more than 2
   * entries in outgoing_particles.
   */
  class InvalidDecay : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };

 private:
  /**
   * Optional success flag for sampling outgoing particles.
   *
   * Set to `true` if 2-body phase-space sampling succeeded, `false` if
   * `sample_masses()` signaled failure via NaNs (e.g., because the
   * chosen channel turned out to be kinematically forbidden after
   * considering potentials), or `std::nullopt` if not sampled yet.
   *
   * This is a **temporary workaround**: the decay channel is currently
   * chosen without knowledge of the potentials, and kinematic failure
   * is handled a posteriori. In the future, this should be replaced by
   * a proper channel selection that already accounts for potential
   * effects during the decision, avoiding the need for this flag.
   */
  std::optional<bool> was_2body_phase_space_sampled_with_potentials_as_valid_ =
      std::nullopt;

 protected:
  /**
   * \ingroup logging
   * Writes information about this decay action to the \p out stream.
   */
  void format_debug_output(std::ostream &out) const override;

  /**
   * Sample outgoing particle types, masses and angles
   * return success of sampling
   */
  bool sample_outgoing_particles();

  /// List of possible decays
  DecayBranchList decay_channels_;

  /// total decay width
  double total_width_;

  /// partial decay width to the chosen outgoing channel
  double partial_width_;

  /// Angular momentum of the decay
  int L_ = 0;

 private:
  /// Spin interaction type
  SpinInteractionType spin_interaction_type_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_DECAYACTION_H_
