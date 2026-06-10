/*
 *
 *    Copyright (c) 2026
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
 * The pion electromagnetic form factor (PEFF) from \iref{Shyam:2010vr} is
 * optionally applied, modifying the differential cross section by
 * \f$|F_\pi(M^2)|^2\f$. The form factor accounts for the internal charged
 * pion propagator.
 *
 * Kinematic variables sampled for the dilepton pair are:
 *   - m_inv:        invariant mass of the dilepton pair
 *                   [\f$2m_e,\sqrt{s}-2m_N\f$]
 *   - q:            3-momentum of dilepton in pn-CM frame
 *   - \f$\cos\theta\f$: Cosine of polar angle of dilepton in pn-CM frame
 *                      \f$[-1, 1]\f$
 *   - \f$\phi\f$:   azimuthal angle of dilepton in pn-CM frame \f$[0, 2\pi]\f$
 *
 * The dilepton 4-momentum is constructed directly from
 * (m_inv, q, \f$\cos\theta\f$, \f$\phi\f$) to enable event-by-event acceptance
 * cuts. The e⁺e⁻ pair is subsequently produced isotropically in the virtual
 * photon's rest frame.
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
  BremsstrahlungActionDilepton(const ParticleList &in, double time,
                               double hadronic_cross_section_input,
                               DileptonBremsPionFormFactor ff_type);

  /**
   * Check if particles can undergo an implemented dilepton
   * bremsstrahlung process.
   *
   * This function does not check the involved kinematics.
   *
   * \param[in] in ParticleList of incoming particles.
   *
   * \return bool if dilepton bremsstrahlung reaction implemented.
   */
  static bool is_dilepton_brems_reaction(const ParticleList &in) {
    return dilepton_brems_reaction_type_(in) != ReactionType::no_reaction;
  }

  /**
   * Adds the hadronic process with a given cross section.
   *
   * The intended use is to add the hadronic cross section from the already
   * performed hadronic action without recomputing it.
   *
   * \param[in] reaction_cross_section Total cross section of underlying
   *                                   hadronic process [mb]
   */
  void add_dummy_hadronic_process(double reaction_cross_section);

  /**
   * Create the final state and write to output.
   *
   * \param[in] outputs List of all outputs. Does not have to be a specific
   *                    dilepton output, the function will take care of this.
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
   *
   * \return ReactionType enum-member
   */
  static ReactionType dilepton_brems_reaction_type_(const ParticleList &in);

  /**
   * Generates momenta of outgoing particles (for 2-body isotropic decays only).
   *
   * \param[in] p_parent FourVector of incoming particle momentum.
   * \param[in] child_1 First child particle.
   * \param[in] child_2 Second child particle.
   *
   * \return    bool if sampling is successful, i.e. if there is enough energy
   *            to create the outgoing particles. If sampling fails, the
   *            function returns false and the weight of the dilepton pair
   *            is set to 0 when generating the final state.
   */
  bool sample_2body_isotropic_(const FourVector &p_parent,
                               ParticleData &child_1, ParticleData &child_2);

  /**
   * Fully differential cross section \f$\frac{d\sigma}{dM dq d\Omega}\f$
   * for \f$ pn \rightarrow pne^+ e^- \f$.
   *
   * \param[in] m_inv  Invariant mass of dilepton pair
   * \param[in] q      3-momentum of dilepton in pn-CM frame
   * \param[in] sqrts  CM energy \f$\sqrt(s)\f$
   *
   * \return           \f$\frac{d\sigma}{dM dq d\Omega}\f$
   */
  double diff_xs_pn_dilepton_(double m_inv, double q, double sqrts) const;

  /**
   * Returns \f$|F_\pi(m_{inv}^2)|^2\f$, the squared pion electromagnetic form
   * factor.
   *
   * \param[in] m_inv_sqr  Invariant mass squared of dilepton pair [GeV²]
   *
   * \return \f$|F_\pi(m_{inv}^2)|^2\f$ (dimensionless)
   */
  double pion_em_form_factor_sq_(double m_inv_sqr) const;

  /**
   * Energy-dependent \f$\Gamma_\rho(m_{inv}^2)\f$ used in the form factor.
   *
   * \param[in] m_inv_sqr  Invariant mass squared [GeV²]
   *
   * \return \f$\Gamma_\rho(m_{inv}^2)\f$ [GeV]
   */
  double gamma_rho_(double m_inv_sqr) const;

  /**
   * Holds the bremsstrahlung branch. As of now, this will hold only one branch.
   */
  CollisionBranchList collision_processes_dilepton_bremsstrahlung_;

  /// Reaction process as determined from incoming particles.
  const ReactionType reaction_type_;

  /// Weight of the dilepton event.
  double weight_ = 0.0;

  /// Total cross section of dilepton bremsstrahlung process [mb].
  double cross_section_dilepton_bremsstrahlung_ = 0.0;

  /// Total hadronic cross section
  const double hadronic_cross_section_;

  /// Form factor type: Off, FF1, FF2.
  const DileptonBremsPionFormFactor form_factor_type_;

  /// Sampled 3-momentum of dilepton pair in pn-CM frame
  double q_;

  /// Sampled invariant mass of the dilepton pair
  double m_inv_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_BREMSSTRAHLUNGACTIONDILEPTON_H_
