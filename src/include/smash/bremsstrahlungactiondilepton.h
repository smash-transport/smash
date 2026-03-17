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
 *   - M  : invariant mass of the dilepton pair  [2me, M_max]
 *   - q  : 3-momentum of dilepton in pn-CM frame
 *   - θ  : polar angle of dilepton in pn-CM frame  [0, π]
 *
 * The dilepton 4-momentum is constructed directly from (M, q, θ) to enable
 * event-by-event acceptance cuts (e.g. HADES filter).
 * The e⁺e⁻ pair is subsequently produced isotropically in the virtual photon's
 * rest frame.
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
   *   3. theta -> polar angle, sampled uniformly in [0, π]
   *   4. phi   -> azimuthal angle, sampled uniformly in [0, 2π]
   *
   * Phasespace construction (two-stage):
   *   Stage 1: Dilepton 4-momentum p_ll built from (M, q, theta, phi) directly.
   *            Nucleon recoil checked.
   *   Stage 2: virtual photon(M) -> e⁺e⁻ isotropically in photon rest frame,
   *            then boosted back to CM frame.
   *
   * Weight set to:
   *   w = dσ/(dM dq dΩ) × ΔM × Δq × 4π² / (σ_had)
   */
  void generate_final_state() override;

  /**
   * Return the weight of the last created photon.
   *
   * \return The total weight.
   */
  double get_total_weight() const override { return weight_; }

  /**
   * Return the total cross section of the underlying hadronic scattering
   * It is necessary for the weighting procedure.
   *
   * \return total cross-section [mb]
   */
  double hadronic_cross_section() const { return hadronic_cross_section_; }

 private:
  // ── Core sampling and kinematics ──────────────────────────────────────────

  /**
   * Generates momenta of outgoing particles (for 2-body isotropic decays only).
   * 
   * \param[in] p_parent FourVector of incoming particle momentum.
   * \param[in] daughter1 First daughter particle.
   * \param[in] daughter2 Second daughter particle.
   */
  void sample_2body_isotropic(const FourVector &p_parent, ParticleData &daughter1, 
    ParticleData &daughter2);

  /**
   * Holds the bremsstrahlung branch. As of now, this will hold only one branch.
   */
  CollisionBranchList collision_processes_dilepton_bremsstrahlung_;

  /// Reaction process as determined from incoming particles.
  const ReactionType reac_;

  /// Form factor type: FF1, FF2 or no_form_factor.
  const FormFactorType form_factor_type_;

  /// Weight of the produced photon.
  double weight_ = 0.0;

  /// Total cross section of dilepton bremsstrahlung process [mb].
  double cross_section_dilepton_bremsstrahlung_ = 0.0;

  /// Total hadronic cross section
  const double hadronic_cross_section_;

  /// Sampled 3-momentum of dilepton pair in pn-CM frame (virtual photon momentum)
  double q_;

  /// Sampled value of theta dilepton in pn-CM frame (polar angle of virtual photon)
  double theta_;

  /// Sampled invariant mass of the dilepton pair
  double M_;

/**
   * Computes the differential cross sections dSigma/dk and dSigma/dtheta of the
   * bremsstrahlung process.
   *
   * \returns Pair containing dSigma/dk as a first argument and dSigma/dtheta
   *          as a second argument
   */
  std::pair<double, double> dilepton_brems_diff_cross_sections();
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_BREMSSTRAHLUNGACTIONDILEPTON_H_
