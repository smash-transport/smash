/*
 *
 *    Copyright (c) 2019-2022,2026
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_BREMSSTRAHLUNGACTIONDILEPTON_H_
#define SRC_INCLUDE_SMASH_BREMSSTRAHLUNGACTIONDILEPTON_H_

#include <utility>

#include "forwarddeclarations.h"
#include "scatteraction.h"

namespace smash {
/**
 * \ingroup action
 * Similar to the photon treatment, BremsstrahlungActionDilepton is a special
 * action which takes two incoming particles and performs a perturbative
 * scattering where a Bremsstrahlung photon is produced.
 * The final state particles are not further propagated, only written
 * to the dilepton output.
 *
 * Implements dilepton production via the process pn -> pn e⁺e⁻ using
 * the phase-space corrected soft-photon approximation (SPA) as e.g. outlined
 * in \iref{Weil:2013mya} (see eq. (40)–(43) in the reference document).
 *
 * The electromagnetic pion form factor (PEFF) from \iref{Shyam:2010vr} is
 * optionally applied, modifying the differential cross section by |F_pi(M²)|².
 * The form factor accounts for the internal charged pion propagator.
 *
 * Kinematic variables sampled (analogous to BremsstrahlungAction):
 *   - M:     invariant mass of the dilepton pair  [2me, M_max]
 *   - q:     3-momentum of dilepton in pn-CM frame
 *   - θ:     polar angle of dilepton in pn-CM frame  [0, pi]
 *   - phi:   azimuthal angle of dilepton in pn-CM frame  [0, 2pi]
 *
 * The dilepton 4-momentum is constructed directly from (M, q, θ, phi) to enable
 * event-by-event acceptance cuts.
 * The e⁺e⁻ pair is subsequently produced isotropically in the virtual photon's
 * rest frame.
 *
 */
class BremsstrahlungActionDilepton : public ScatterAction {
 public:
  /**
   * Enum for encoding bremsstrahlung process for n+p only.
   * It is uniquely determined by the incoming particles.
   */
  enum class ReactionType { no_reaction, np };

  /**
   * Construct a BremsstrahlungActionDilepton object.
   *
   * \param[in] in ParticleList of incoming particles (n+p only).
   * \param[in] time Time relative to underlying hadronic action.
   * \param[in] hadronic_cross_section_input Total np hadronic cross section.
   * \param[in] ff_type Which form factor parametrization to use.
   *
   * \return The constructed object.
   */
  BremsstrahlungActionDilepton(const ParticleList &in, const double time,
                               const double hadronic_cross_section_input,
                               DileptonBremsPionFormFactor ff_type);

  /**
   * Check if particles can undergo an implemented dilepton
   * bremsstrahlung process for invocation in experiment.h/.cc.
   *
   * This function does not check the involved kinematics.
   *
   * \param[in] in ParticleList of incoming particles.
   * \return bool if dilepton bremsstrahlung reaction implemented.
   */
  static bool is_dilepton_brems_reaction(const ParticleList &in) {
    return dilepton_brems_reaction_type(in) != ReactionType::no_reaction;
  }

  /**
   * Adds the hadronic process with a given cross section.
   *
   * The intended use is to add the hadronic cross section from the already
   * performed hadronic action without recomputing it.
   *
   * \param[in] reaction_cross_section Total cross section of underlying
   *                                    hadronic process [mb]
   */
  void add_dummy_hadronic_process(double reaction_cross_section);

  /**
   * Create the final state and write to output.
   *
   * \param[in] outputs List of all outputs. Does not have to be a specific
   *                      dilepton output, the function will take care of this.
   */
  void perform_dilepton_bremsstrahlung(const OutputsList &outputs);

  /**
   * Main function: sample kinematics and compute weight for one
   * dilepton event.
   */
  void generate_final_state() override;

  /**
   * Return the weight of the dilepton pair.
   *
   * \return The total weight.
   */
  double get_total_weight() const override { return weight_; }

  /**
   * Return the partial weight of zero as otherwise garbage is written to
   * output.
   *
   * \return 0.0 as partial weight.
   */
  double get_partial_weight() const override { return 0.0; }

 private:
  /**
   * Determine dilepton bremsstrahlung process from incoming particle list.
   *
   * If incoming particles are not part of any implemented dilepton
   * bremsstrahlung process, i.e. currently only p+n (in any order),
   * the function will return no_reaction.
   *
   * \param[in] in ParticleList of incoming particles.
   * \return ReactionType enum-member
   */
  static ReactionType dilepton_brems_reaction_type(const ParticleList &in);

  /**
   * Generates momenta of outgoing particles (for 2-body isotropic decays only).
   *
   * \param[in] p_parent FourVector of incoming particle momentum.
   * \param[in] daughter1 First daughter particle.
   * \param[in] daughter2 Second daughter particle.
   * \return    bool if sampling successful, i.e. if there is enough energy
   *            to create the outgoing particles. If sampling fails, the
   * function returns false and the weight of the dilepton pair is set to 0 in
   * generate_final_state().
   */
  bool sample_2body_isotropic(const FourVector &p_parent,
                              ParticleData &daughter1, ParticleData &daughter2);

  /**
   * Fully differential cross section dsigma/(dM dq dOmega) for pn -> pne⁺e⁻.
   *
   * \param[in] M      Invariant mass of dilepton pair
   * \param[in] q      3-momentum of dilepton in pn-CM frame
   * \param[in] sqrts  CM energy sqrt(s)
   * \return           dsigma/(dM dq dOmega)
   */
  double diff_xs_pn_dilepton(double M, double q, double sqrts) const;

  /**
   * Returns |F_pi(M²)|², the squared pion electromagnetic form factor.
   *
   * \param[in] M2  M² = invariant mass squared of dilepton pair [GeV²]
   * \return |F_pi(M²)|² (dimensionless)
   */
  double pion_em_form_factor_sq(double M2) const;

  /**
   * Energy-dependent Gamma_rho(M²) used in the form factor.
   *
   * \param[in] M_sq  Invariant mass squared [GeV²]
   * \return Gamma_rho(M²) [GeV]
   */
  double gamma_rho(double M_sq) const;

  /**
   * Helper function for calculating R_2 as defined in \iref{Weil:2013mya},
   * eq. (42).
   *
   * \param[in] s Mandelstam variable s [GeV²]
   * \return R_2(M²) (dimensionless)
   */
  double R_2_helper(const double s) const;
  /**
   * Holds the bremsstrahlung branch. As of now, this will hold only one branch.
   */
  CollisionBranchList collision_processes_dilepton_bremsstrahlung_;

  /// Reaction process as determined from incoming particles.
  const ReactionType reac_;

  /// Weight of the dilepton event.
  double weight_ = 0.0;

  /// Total cross section of dilepton bremsstrahlung process [mb].
  double cross_section_dilepton_bremsstrahlung_ = 0.0;

  /// Total hadronic cross section
  const double hadronic_cross_section_;

  /// Form factor type: Off (FF=1), FF1, FF2.
  const DileptonBremsPionFormFactor form_factor_type_;

  /// Sampled 3-momentum of dilepton pair in pn-CM frame
  double q_;

  /// Sampled value of dilepton polar angle in pn-CM frame
  double theta_;

  /// Sampled invariant mass of the dilepton pair
  double M_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_BREMSSTRAHLUNGACTIONDILEPTON_H_
