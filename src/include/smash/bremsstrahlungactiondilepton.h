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
  // ──────── Constructor ─────────────────────────── 
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
  // ──────── Further nonstatic methods ───────────────────────────────
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
   * Compute dilepton bremsstrahlung cross sections and register
   * the single pn -> pne⁺e⁻ process branch. Must be called after 
   * add_dummy_hadronic_process() and before perform_dilepton_bremsstrahlung().
   * Analogous to the photon bremsstrahlung add_single_process().
   */
  void add_single_process() {
    add_processes<CollisionBranch>(dilepton_brems_cross_sections(),
                                   collision_processes_dilepton_bremsstrahlung_,
                                   cross_section_dilepton_bremsstrahlung_);
  }

  /**
   * Create the final state and write to output.
   *
   * \param[in] outputs List of all outputs. Does not have to be a specific
   *                      dilepton output, the function will take care of this.
   */
  void perform_dilepton_bremsstrahlung(const OutputsList &outputs);

  /**
   * Generate the final-state for the Bremsstrahlung process. Generates only
   * 3-body final state.
   */
  void generate_final_state() override;

  /**
   * Sample the final state anisotropically, considering the differential
   * cross sections with respect to theta and k.
   */
  void sample_3body_phasespace();

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
  /**
   * Holds the bremsstrahlung branch. As of now, this will hold only one branch.
   */
  CollisionBranchList collision_processes_dilepton_bremsstrahlung_;

  /// Reaction process as determined from incoming particles.
  const ReactionType reac_;

  /// Form factor type: FF1, FF2 or no_form_factor.
  const FormFactorType form_factor_type_;

  //TODO: Probably not needed anymore but to be checked after implementation.
  /// Weight of the produced photon.
  /// double weight_ = 0.0;

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
   * Computes the total cross section of the dilepton bremsstrahlung process.
   *
   * \returns List of photon reaction branches.
   */
  CollisionBranchList dilepton_brems_cross_sections();

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
