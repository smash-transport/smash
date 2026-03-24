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

#include "scatteraction.h"
#include "forwarddeclarations.h"

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
 * in Weil (2013) (see Eq. (40)–(43) in the reference document).
 *
 * The electromagnetic pion form factor (PEFF) from Shyam & Mosel (2010),
 * Phys. Rev. C 82, 062201(R) is optionally applied, modifying the
 * differential cross section by |F_\pi(M²)|². The form factor accounts
 * for the internal charged pion propagator in diagram 1(c) (pn-channel
 * only, charged meson exchange).
 *
 * Kinematic variables sampled (analogous to BremsstrahlungAction):
 *   - M:     invariant mass of the dilepton pair  [2me, M_max]
 *   - q:     3-momentum of dilepton in pn-CM frame
 *   - theta: polar angle of dilepton in pn-CM frame  [0, pi]
 *   - phi:   azimuthal angle of dilepton in pn-CM frame  [0, 2pi]
 *
 * The dilepton 4-momentum is constructed directly from (M, q, θ) to enable
 * event-by-event acceptance cuts (e.g. HADES filter).
 * The e⁺e⁻ pair is subsequently produced isotropically in the virtual photon's
 * rest frame.
 * 
 * Inherited member from Action (3rd level: Action -> ScatterAction -> here):
 * weight_              : Weight of the produced dilepton pair
 * outgoing_particles_  : Final state particles {p, n, e⁺, e⁻}
 * time_of_execution_   : Time of execution
 * incoming_particles_  : Incoming particles {p, n}
 */
class BremsstrahlungActionDilepton : public ScatterAction {
 public:
  /**
   * Enum for encoding bremsstrahlung process for n+p only. 
   * It is uniquely determined by the incoming particles. 
   * The naming scheme is :
   * Incoming_1_Incoming_2_.
   */
  enum class ReactionType {no_reaction, np};

  /**
   * Choice of pion electromagnetic form factor parametrization acc. to
   * Shyam & Mosel (2010).:
   * FF1: photon couples to pion via pure \rho0 meson
   * FF2: photon couples 50% directly to intrinsic quark structure of pion
   *      and 50% indirectly via \rho0 meson
   */
  enum class FormFactorType {FF1, FF2, no_form_factor};

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
  // ──────── Constructor ─────────────────────────────────────────────────────── 
  BremsstrahlungActionDilepton(const ParticleList &in, const double time,
                       const double hadronic_cross_section_input, 
                       FormFactorType ff_type);
  
  // ──────── Static methods for preprocessing checks ─────────────────────────── 
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
  // ──────── Further nonstatic methods ─────────────────────────────────────────
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
   *
   * Sampling strategy (four variables, analogous to k and theta in photons):
   *   1. M     -> invariant mass, sampled uniformly in [2*m_e, M_max]
   *   2. q     -> 3-momentum, sampled uniformly in [q_min, q_max], limits depend 
   *               on M and sqrt_s
   *   3. theta -> polar angle, sampled uniformly in [0, pi]
   *   4. phi   -> azimuthal angle, sampled uniformly in [0, 2*pi]
   *
   * Phasespace construction (two-stage):
   *   Stage 1: Dilepton 4-momentum p_ll built from (M, q, theta, phi) directly.
   *            Nucleon recoil checked.
   *   Stage 2: virtual photon(M) -> e⁺e⁻ isotropically in photon rest frame,
   *            then boosted back to CM frame.
   *
   * Weight set to:
   *   w = dsigma/(dM dq dOmega) * Delta_M * Delta_q * 2*pi² / (sigma_hadronic)
   */
  void generate_final_state() override;

  // ── Add override function to handle weight retrieval ────────────────────────
  /**
   * Return the weight of the last created photon.
   *
   * \return The total weight.
   */
  double get_total_weight() const override { return weight_; }

  // ── Add override function to handle partial weight ──────────────────────────
  /**
   * Return the partial weight of zero as otherwise garbage was written to output.
   *
   * \return 0.0 as partial weight.
   */
  double get_partial_weight() const override { return 0.0; }

 private:
  // ── Core sampling and kinematics ────────────────────────────────────────────
  /**
   * Generates momenta of outgoing particles (for 2-body isotropic decays only).
   * 
   * \param[in] p_parent FourVector of incoming particle momentum.
   * \param[in] daughter1 First daughter particle.
   * \param[in] daughter2 Second daughter particle.
   * \return    bool if sampling successful, i.e. if there is enough energy 
   *            to create the outgoing particles. If sampling fails, the function 
   *            returns false and the weight of the dilepton pair is set to 0 
   *            in generate_final_state().
   *            Note: This function is used for the "decay" of the virtual photon 
   *                  into e⁺e⁻, but can be used for any isotropic 2-body decay.
   */
  bool sample_2body_isotropic(const FourVector &p_parent, ParticleData &daughter1, 
    ParticleData &daughter2);

  // ── Differential cross section: SPA + PEFF ──────────────────────────────────
  /**
   * Fully differential cross section dsigma/(dM dq dtheta) for pn -> pne⁺e⁻.
   *
   * Based on the phase-space corrected SPA:
   * dsigma/(dM dE dOmega) = (alpha²/6pi³) * (q/ME²) * sigma(s) * R2(s2)/R2(s)
   * with:
   *   sigma(s)  = [s - (m_p + m_n)²] / (2*m_p²) * sigma_el(s)
   *   R2(s) = sqrt(1 - (m_p + m_n)²/s)
   *   s2    = s + M² - 2*E*sqrt(s)
   * 
   * Variable substitution: E = sqrt(q² + M²), so dE = q/E dq,
   * resulting in dsigma/(dM dq dOmega) = dsigma/(dM dE dOmega) * (q/E).
   *
   * \param[in] M      Invariant mass of dilepton pair
   * \param[in] q      3-momentum of dilepton in pn-CM frame
   * \param[in] sqrts  CM energy sqrt(s)
   * \return           dsigma/(dM dq dOmega)
   */
  double diff_xs_pn_dilepton(double M, double q, double sqrts) const;

  // ── Pion electromagnetic form factor (Shyam & Mosel 2010) ───────────────────
  /**
   * Returns |F_pi(M²)|², the squared pion electromagnetic form factor.
   *
   * The energy-dependent ρ width Γ_ρ(M²) follows the standard
   * parametrization already available in SMASH (decaytype.cc).
   *
   * \param[in] M2  M² = invariant mass squared of dilepton pair [GeV²]
   * \return |F_pi(M²)|² (dimensionless)
   */
  double pion_em_form_factor_sq(double M2) const;

  // ── Gamma_rho (Shyam & Mosel 2010) ──────────────────────────────────────────
  /**
   * Energy-dependent Gamma_rho(M²) used in the form factor.
   *
   * \param[in] M_sq  Invariant mass squared [GeV²]
   * \return Gamma_rho(M²) [GeV]
   */
  double gamma_rho(double M_sq) const;

  // ── Little helper function ──────────────────────────────────────────────────
  double R_2_helper(const double s) const;
  /**
   * Holds the bremsstrahlung branch. As of now, this will hold only one branch.
   */
  CollisionBranchList collision_processes_dilepton_bremsstrahlung_;

  /// Reaction process as determined from incoming particles.
  const ReactionType reac_;

  /// Weight of the dilepton event.
  /// Set in generate_final_state(), used by perform_dilepton_bremsstrahlung()
  /// for the output. Analogous to weight_ in BremsstrahlungAction. Requires
  /// the override of get_total_weight() to work with the output machinery.
  double weight_ = 0.0;

  /// Total cross section of dilepton bremsstrahlung process [mb].
  double cross_section_dilepton_bremsstrahlung_ = 0.0;

  /// Total hadronic cross section
  const double hadronic_cross_section_;

  /// Form factor type: FF1, FF2 or no_form_factor.
  const FormFactorType form_factor_type_;

  /// Sampled 3-momentum of dilepton pair in pn-CM frame (virtual photon momentum)
  double q_;

  /// Sampled value of theta dilepton in pn-CM frame (polar angle of virtual photon)
  double theta_;

  /// Sampled invariant mass of the dilepton pair
  double M_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_BREMSSTRAHLUNGACTIONDILEPTON_H_
